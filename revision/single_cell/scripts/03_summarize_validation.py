"""Patient-level validation statistics; plotting scripts only read these tables."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon, rankdata
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'data/processed'
RNG=np.random.default_rng(20260905)
V=pd.read_csv(P/'validation_primary_scores.csv')
D=pd.read_csv(P/'discovery_validation_ready.csv')
PROGRAMS=['MP1','MP2','MP3','MP4','MP4_without_ROS']

def boot_mean_ci(x,B=5000):
    x=np.asarray(x,dtype=float)
    if len(x)<3:return np.nan,np.nan
    b=x[RNG.integers(0,len(x),(B,len(x)))].mean(axis=1)
    return tuple(np.quantile(b,[.025,.975]))

def signed_p(x):
    x=np.asarray(x,dtype=float)
    if len(x)<3:return np.nan
    if np.all(x==0):return 1.
    return float(wilcoxon(x,alternative='two-sided',method='auto').pvalue)

all_pairs=[];tests=[];corr=[];technical=[]
for dataset,df,idcol,unit in [('GSE149614',V,'patient','patient'),('GSE149614_QC_filtered',V[V.technical_qc_pass.astype(bool)].copy(),'patient','patient'),('GSE202642',D,'sample','tissue sample')]:
    for gid,g in df.groupby(idcol):
        for p in PROGRAMS:
            score=g[p+'_UCell']
            rho=float(spearmanr(g.OS_composite,score).statistic)
            # Partial rank correlation controls technical covariates within sample.
            mt=g['percent.mt'] if 'percent.mt' in g else pd.Series(np.zeros(len(g)),index=g.index)
            X=np.column_stack([np.ones(len(g)),rankdata(np.log1p(g.nCount_RNA)),rankdata(mt)])
            xr=rankdata(g.OS_composite);yr=rankdata(score)
            xr=xr-X@np.linalg.lstsq(X,xr,rcond=None)[0];yr=yr-X@np.linalg.lstsq(X,yr,rcond=None)[0]
            partial=float(np.corrcoef(xr,yr)[0,1]) if np.std(yr)>0 else np.nan
            corr.append(dict(dataset=dataset,unit=unit,sample=gid,program=p,n_cells=len(g),rho=rho,partial_rho_depth_mt=partial))
    definitions=['OS_high_top20','OS_high_within_sample_top20']
    if 'OS_UCell_high_within_sample_top20' in df:definitions+=['OS_UCell_high_within_sample_top20']
    for definition in definitions:
        for prog in PROGRAMS:
            pairs=[]
            for gid,g in df.groupby(idcol):
                low=g.loc[~g[definition].astype(bool),prog+'_UCell'];high=g.loc[g[definition].astype(bool),prog+'_UCell']
                eligible=len(low)>=5 and len(high)>=5
                row=dict(dataset=dataset,unit=unit,definition=definition,program=prog,sample=gid,n_low=len(low),n_high=len(high),mean_low=low.mean(),mean_high=high.mean(),difference=high.mean()-low.mean(),eligible=eligible)
                all_pairs.append(row)
                if eligible:pairs.append(row['difference'])
            ci=boot_mean_ci(pairs)
            tests.append(dict(dataset=dataset,unit=unit,definition=definition,program=prog,n_paired=len(pairs),mean_difference=np.mean(pairs) if pairs else np.nan,ci_low=ci[0],ci_high=ci[1],P=signed_p(pairs)))
    technical.append(dict(dataset=dataset,n_cells=len(df),n_units=df[idcol].nunique(),unit=unit))

t=pd.DataFrame(tests)
t['FDR_within_definition']=np.nan
for _,ix in t.groupby(['dataset','definition']).groups.items():
    valid=t.loc[ix,'P'].dropna()
    if len(valid):t.loc[valid.index,'FDR_within_definition']=multipletests(valid,method='fdr_bh')[1]
t.to_csv(P/'patient_level_program_tests.csv',index=False)
pd.DataFrame(all_pairs).to_csv(P/'patient_level_program_pairs.csv',index=False)
c=pd.DataFrame(corr);c.to_csv(P/'within_patient_correlations.csv',index=False)
ct=[]
for (dataset,prog),g in c.groupby(['dataset','program']):
    for metric in ['rho','partial_rho_depth_mt']:
        x=g[metric].dropna().to_numpy();ci=boot_mean_ci(x)
        ct.append(dict(dataset=dataset,program=prog,metric=metric,n_units=len(x),mean_rho=np.mean(x),ci_low=ci[0],ci_high=ci[1],P=signed_p(x)))
pd.DataFrame(ct).to_csv(P/'patient_correlation_summary.csv',index=False)
pd.DataFrame(technical).to_csv(P/'single_cell_cohort_counts.csv',index=False)
ds=D.groupby('sample').agg(n_cells=('cell','size'),n_high=('OS_high_top20','sum')).reset_index()
ds['fraction_high']=ds.n_high/ds.n_cells
ds.to_csv(P/'discovery_sample_group_counts.csv',index=False)
pd.crosstab(D.OS_group_valley,D.OS_high_top20).reindex(index=['OS-low','OS-high'],columns=[False,True],fill_value=0).to_csv(P/'discovery_group_agreement.csv')
s=D.groupby('sample').agg(n_cells=('cell','size'),historical_high=('OS_group_valley',lambda x:(x=='OS-high').sum()),top20_high=('OS_high_top20','sum'),tissue_group=('group','first')).reset_index()
s['historical_fraction']=s.historical_high/s.n_cells
s['top20_fraction']=s.top20_high/s.n_cells
s=s.iloc[s['sample'].str.extract(r'(\d+)')[0].astype(int).argsort()]
s.to_csv(P/'discovery_historical_and_top20_counts.csv',index=False)
edges=np.linspace(D.OS_composite.min(),D.OS_composite.max(),65)
hist_rows=[]
for lab,g in D.groupby('OS_group_valley'):
    counts,_=np.histogram(g.OS_composite,bins=edges)
    for lo,hi,n in zip(edges[:-1],edges[1:],counts):hist_rows.append(dict(group=lab,bin_left=lo,bin_right=hi,n_cells=n))
pd.DataFrame(hist_rows).to_csv(P/'discovery_score_histogram.csv',index=False)

# Explicit score/cutoff/version record: no silent replacement of historical labels.
audit=pd.read_csv(P/'discovery_definition_audit.csv')
report={'bootstrap_unit':'patient in validation; tissue sample in discovery','bootstrap_replicates':5000,'seed':20260905,'minimum_cells_per_group_per_unit':5,'primary_group_definition':'cohort-wide 80th percentile of composite OS score','sensitivity_group_definitions':['OS_high_within_sample_top20','OS_UCell_high_within_sample_top20'],'technical_QC_sensitivity':'retains original scores and group labels, excludes cells with <200 or >6000 detected genes or >=20% mitochondrial counts','primary_program':'frozen discovery MP4','controls':'MP4 excluding ROS-overlap genes; partial rank correlation controlling library size and mitochondrial percentage','historical_discovery_labels':'retained in OS_group_valley; top20 labels are a new sensitivity/revision definition and must not be described as the historical model derivation labels','discovery_definition_audit':audit.to_dict('records')}
(ROOT/'audit/statistical_protocol.json').write_text(json.dumps(report,indent=2))
print(t[(t.dataset=='GSE149614')&t.program.isin(['MP4','MP4_without_ROS'])].to_string(index=False))
print('VALIDATION_STATS_COMPLETE')
