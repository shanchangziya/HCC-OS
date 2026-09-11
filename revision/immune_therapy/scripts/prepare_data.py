"""Create immutable plotting tables from local patient data and frozen OSARS scores.
Statistical calculations belong here, never in plot_panels.py.
"""
from pathlib import Path
import hashlib, json, os, platform, shutil
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
PROJECT = Path(os.environ.get('HCC_OS_LEGACY_PROJECT_ROOT',ROOT.parents[1])).expanduser().resolve()
INPUT = Path(os.environ.get('HCC_OS_IMMUNE_INPUT_DIR',ROOT/'data/raw')).expanduser().resolve()
OUT = ROOT / 'data/processed'
SOURCE = ROOT / 'data/source'
SEED = 20260905
RNG = np.random.default_rng(SEED)
SOURCES = []

def read(path, **kwargs):
    path = Path(path)
    try:path_label=str(path.resolve().relative_to(ROOT.parents[1]))
    except ValueError:path_label=f'external:{path.name}'
    SOURCES.append({'path':path_label,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size})
    return pd.read_csv(path, **kwargs)

def save(df, name):
    df.to_csv(OUT/name,index=False)

def star(p):
    return '***' if p < .001 else '**' if p < .01 else '*' if p < .05 else 'ns'

def bh(p):
    return stats.false_discovery_control(np.array(p),method='bh')

def stats_cont(long, family, bootstrap=True):
    rows=[]; boxes=[]
    for feature, d in long.groupby('Feature',sort=False):
        lo = d.loc[d.Risk_Group=='Low','Value'].dropna().to_numpy()
        hi = d.loc[d.Risk_Group=='High','Value'].dropna().to_numpy()
        pooled_sd = np.std(np.r_[lo,hi],ddof=1)
        eff=(hi.mean()-lo.mean())/pooled_sd if pooled_sd>0 else np.nan
        ci=[np.nan,np.nan]
        if bootstrap and pooled_sd>0:
            samples=(RNG.choice(hi,(5000,len(hi))).mean(axis=1)-RNG.choice(lo,(5000,len(lo))).mean(axis=1))/pooled_sd
            ci=np.quantile(samples,[.025,.975])
        p=stats.mannwhitneyu(hi,lo,alternative='two-sided',method='asymptotic',use_continuity=True).pvalue if pooled_sd>0 else 1.0
        rows.append(dict(Feature=feature,Family=family,Low_N=len(lo),High_N=len(hi),Low_mean=lo.mean(),High_mean=hi.mean(),
                         Low_median=np.median(lo),High_median=np.median(hi),Difference=hi.mean()-lo.mean(),P_value=p,
                         Standardized_difference=eff,CI_low=ci[0],CI_high=ci[1],Pooled_SD=pooled_sd,
                         Difference_CI_low=ci[0]*pooled_sd,Difference_CI_high=ci[1]*pooled_sd,
                         Test='Two-sided Mann–Whitney U; asymptotic; continuity correction; ties corrected'))
        for group,v in [('Low',lo),('High',hi)]:
            q1,med,q3=np.quantile(v,[.25,.5,.75]);iqr=q3-q1
            wl=v[v>=q1-1.5*iqr].min();wh=v[v<=q3+1.5*iqr].max()
            boxes.append(dict(Feature=feature,Risk_Group=group,N=len(v),q1=q1,med=med,q3=q3,whislo=wl,whishi=wh,
                              vmin=v.min(),vmax=v.max()))
    s=pd.DataFrame(rows);s['FDR']=bh(s.P_value);s['Significance']=s.FDR.map(star)
    return s,pd.DataFrame(boxes)

def attach_jitter(long):
    long=long.copy();long['Jitter']=RNG.uniform(-.18,.18,len(long))
    return long

def wilson(k,n):
    z=stats.norm.ppf(.975);p=k/n;den=1+z*z/n
    m=(p+z*z/(2*n))/den;h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return m-h,m+h

def response_stats(d,cohort):
    rows=[]
    counts=[]
    for g in ['Low','High']:
        dg=d[d.Risk_Group==g]; n=len(dg);k=int((dg.Response=='R').sum());ci=wilson(k,n)
        rows.append(dict(Cohort=cohort,Risk_Group=g,N=n,R=k,NR=n-k,Proportion=k/n,CI_low=ci[0],CI_high=ci[1]))
        counts.append([k,n-k])
    odds,p=stats.fisher_exact([counts[1],counts[0]],alternative='two-sided')
    for r in rows:r.update(P_value=p,Odds_ratio_High_vs_Low=odds,Significance=star(p),Test='Two-sided Fisher exact',CI_type='Wilson 95%')
    return pd.DataFrame(rows)

OUT.mkdir(parents=True,exist_ok=True)
groups=read(PROJECT/'FIG5_redone_drug_mutation/TCGA_LIHC_OSARS_groups.csv')
assert len(groups)==343 and groups.Sample.is_unique
median=groups.Risk_Score.median()
assert (groups.Risk_Group==np.where(groups.Risk_Score>=median,'High','Low')).all()
assert groups.Risk_Group.value_counts().to_dict()=={'High':172,'Low':171}
save(groups,'tcga_frozen_groups.csv')
strict=groups.copy();strict['Risk_Group']=np.where(strict.Risk_Score>median,'High','Low')
boundary=groups.merge(strict,on=['Sample','Risk_Score'],suffixes=('_frozen','_strict'))
save(boundary[boundary.Risk_Group_frozen!=boundary.Risk_Group_strict],'median_boundary_patient.csv')

ss=read(SOURCE/'ssgsea_patient_export.csv')
cib=read(INPUT/'CIBERSORT_Results.txt',sep='\t');cib=cib.rename(columns={cib.columns[0]:'Sample'})
est=read(INPUT/'ESTIMATE_score.txt',sep='\t');est=est.rename(columns={est.columns[0]:'Sample'})
feature_sets=[(ss,{'Activated CD8 T cell':'Activated CD8 T cells','Regulatory T cell':'Regulatory T cells',
                   'Natural killer cell':'Natural killer cells','Macrophage':'Macrophage signature','MDSC':'MDSC signature'},'ssGSEA'),
              (cib,{'T cells regulatory (Tregs)':'Tregs','Macrophages M0':'Macrophages M0','Macrophages M2':'Macrophages M2'},'CIBERSORT'),
              (est,{'ImmuneScore':'Immune score','StromalScore':'Stromal score','TumorPurity':'Tumor purity'},'ESTIMATE')]
immune=[]
for frame,mapping,method in feature_sets:
    for key,label in mapping.items():
        d=groups.merge(frame[['Sample',key]],on='Sample',validate='one_to_one').rename(columns={key:'Value'})
        d['Feature']=label;d['Method']=method;d['Source_feature']=key;immune.append(d)
immune=attach_jitter(pd.concat(immune,ignore_index=True))
assert immune.Value.notna().all() and len(immune)==343*11
ist,ibox=stats_cont(immune,'11 selected immune/stromal features')
ist=ist.merge(immune[['Feature','Method','Source_feature']].drop_duplicates(),on='Feature')
save(immune,'immune_patient_values.csv');save(ist,'immune_statistics.csv');save(ibox,'immune_box_summary.csv')

tide=read(INPUT/'TIDE_TCGA343.txt',sep='\t');tide=tide.rename(columns={tide.columns[0]:'Sample'})
tide=groups.merge(tide,on='Sample',validate='one_to_one');assert len(tide)==343
features=['TIDE','Dysfunction','Exclusion','MDSC','CAF','TAM M2']
tl=attach_jitter(tide.melt(id_vars=['Sample','Risk_Group','Risk_Score'],value_vars=features,var_name='Feature',value_name='Value'))
ts,tb=stats_cont(tl,'Six TIDE-derived scores')
save(tl,'tide_patient_values.csv');save(ts,'tide_statistics.csv');save(tb,'tide_box_summary.csv')
tide['Response']=np.where(tide.Responder.astype(str).str.lower().isin(['true','1','responder']),'R','NR')
save(tide[['Sample','Risk_Group','Risk_Score','Response']],'tide_predicted_response_patients.csv')
tresp=response_stats(tide,'TCGA-LIHC: TIDE predicted ICB response')
assert list(tresp.R)==[84,55] and list(tresp.N)==[171,172]
save(tresp,'tide_response_statistics.csv')

# ACT is repredicted from the frozen GBM with canonical preprocessing.
actaudit=read(SOURCE/'ACT_input_audit.csv').iloc[0]
assert actaudit.Duplicate_normalized_genes==0 and actaudit.NA_cells==0 and actaudit.Missing_genes=='GC'
legacy=read(PROJECT/'FIG5_therapy_immune/Fig5F_external_treatment_OSARS_scores.csv')
act=read(SOURCE/'ACT_scores_current_model.csv')
old_act=legacy[legacy.Cohort.str.contains('GSE100797')]
score_check=act.merge(old_act[['ID','Risk_Score']],on='ID',suffixes=('_new','_old'))
assert np.allclose(score_check.Risk_Score_new,score_check.Risk_Score_old,atol=1e-10,rtol=0)
clin=read(SOURCE/'ACT_clinical_export.csv')
assert len(act)==21 and (act.RR==np.where(act.Response.isin(['CR','PR']),'R','NR')).all()
act=act.rename(columns={'Response':'RECIST','RR':'Response'})
act['Risk_Group']=np.where(act.Risk_Score>act.Risk_Score.median(),'High','Low')
act['Cohort']='Melanoma TIL-ACT (GSE100797)'
assert act.Risk_Group.value_counts().to_dict()=={'Low':11,'High':10}

therapy=[];therstats=[]
for filename,cohort in [('tace_OSARS_scores_helper_train_mean.csv','TACE (GSE104580)'),
                        ('sorafenib_OSARS_scores_helper_train_mean.csv','Adjuvant sorafenib (GSE109211)')]:
    d=read(PROJECT/'FIG5_therapy_recomputed_20260904'/filename)
    assert d.ID.is_unique and (d.Risk_Group==np.where(d.Risk_Score>d.Risk_Score.median(),'High','Low')).all()
    d['Cohort']=cohort;therapy.append(d);therstats.append(response_stats(d,cohort))
therapy.append(act);therstats.append(response_stats(act,act.Cohort.iloc[0]))
txstats=pd.concat(therstats,ignore_index=True)
q=bh(txstats.drop_duplicates('Cohort').P_value)
txstats['FDR_three_exploratory_cohorts']=txstats.Cohort.map(dict(zip(txstats.Cohort.unique(),q)))
save(pd.concat(therapy,ignore_index=True),'therapy_patients.csv');save(txstats,'therapy_response_statistics.csv')
ref=read(PROJECT/'FIG5_therapy_recomputed_20260904/TACE_sorafenib_recomputed_summary_helper.csv')
for idx,cohort in enumerate(txstats.Cohort.unique()[:2]):
    rs=txstats[txstats.Cohort==cohort]
    assert np.isclose(rs.P_value.iloc[0],ref.Fisher_two_sided_P.iloc[idx],rtol=1e-10)
    assert list(rs.R)==[ref.Low_R.iloc[idx],ref.High_R.iloc[idx]]

# Strict-median sensitivity: recompute groups and all original feature-family FDRs.
sens=[]
for label,long in [('immune',immune),('tide',tl)]:
    sl=long.drop(columns=['Risk_Group','Risk_Score']).merge(strict,on='Sample',validate='many_to_one')
    st,_=stats_cont(sl,label+' strict median sensitivity',bootstrap=False);st['Module']=label;sens.append(st)
save(pd.concat(sens,ignore_index=True),'strict_median_sensitivity.csv')
sp=tide.drop(columns=['Risk_Group','Risk_Score']).merge(strict,on='Sample',validate='one_to_one')
save(response_stats(sp,'TCGA strict median sensitivity'),'strict_median_tide_response.csv')

# CIBERSORT fit sensitivity includes all 22 fractions, with a separate BH family.
cib_sensitivity=[]
for include,frame in [('all_343',cib),('CIBERSORT_P_lt_0.05',cib[cib['P-value']<.05])]:
    cc=groups.merge(frame,on='Sample',validate='one_to_one')
    cl=cc.melt(id_vars=['Sample','Risk_Group'],value_vars=list(cib.columns[1:23]),var_name='Feature',value_name='Value')
    cs,_=stats_cont(cl,include+';22 cell fractions',bootstrap=True);cs['Filter']=include;cib_sensitivity.append(cs)
save(pd.concat(cib_sensitivity,ignore_index=True),'cibersort_22_fraction_sensitivity.csv')

drug=read(PROJECT/'FIG5_integrated_final/DrugReflector_OSARS_High_to_Low_top100.csv')
anno=read(PROJECT/'FIG5_integrated_final/DrugReflector_OSARS_High_to_Low_top18_annotated.csv')
drug=drug.merge(anno[['compound','label']],on='compound',how='left',validate='one_to_one')
drug['Label']=drug.label.fillna(drug.compound)
drug['Model_probability_percent']=drug.probability*100
drug['Rank_1based']=drug['rank'].astype(int)+1
save(drug,'drugreflector_top100.csv');save(drug.sort_values('rank').head(20),'drugreflector_top20.csv')

genes=['PDCD1','CD274','CTLA4','LAG3','TIGIT','HAVCR2','VSIR','IDO1']
expr_path=INPUT/'TCGA343_exprset.txt'
try:expr_label=str(expr_path.resolve().relative_to(ROOT.parents[1]))
except ValueError:expr_label=f'external:{expr_path.name}'
SOURCES.append({'path':expr_label,'sha256':hashlib.sha256(expr_path.read_bytes()).hexdigest(),'bytes':expr_path.stat().st_size})
parts=[]
for chunk in pd.read_csv(expr_path,sep='\t',chunksize=1000):
    parts.append(chunk[chunk.iloc[:,0].str.strip().isin(genes)])
expr=pd.concat(parts);expr=expr.rename(columns={expr.columns[0]:'Feature'})
ck=expr.melt(id_vars='Feature',var_name='Sample',value_name='Value').merge(groups,on='Sample',validate='many_to_one')
ck=attach_jitter(ck);cst,cbox=stats_cont(ck,'Eight checkpoint genes',bootstrap=False)
save(ck,'checkpoint_patient_values.csv');save(cst,'checkpoint_statistics.csv');save(cbox,'checkpoint_box_summary.csv')

manifest={'seed':SEED,'bootstrap_replicates':5000,'tcga_median':median,'frozen_boundary':'High >= cohort median; Low < median',
          'external_boundary':'High > eligible-cohort median; Low <= median','sources':SOURCES,
          'software':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__},
          'ACT_score_source':'Repredicted frozen GBM with canonical helper preprocessing; GC imputed by training mean; > median restores Article2 11/10 grouping; scores agree with archive within 1e-10',
          'clinical_endpoint_correction':'BIOSTORM labels denote the original recurrence-benefit gene-signature class, not radiographic tumor shrinkage; ACT labels are RECIST CR/PR versus SD/PD in melanoma'}
(OUT/'analysis_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print(ist[['Feature','Difference','FDR']].to_string(index=False))
print(ts[['Feature','Difference','FDR']].to_string(index=False))
print(txstats[['Cohort','Risk_Group','R','N','Proportion','P_value','FDR_three_exploratory_cohorts']].to_string(index=False))
print(tresp.to_string(index=False))
