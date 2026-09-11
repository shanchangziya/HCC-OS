"""Generate plotting-ready score, discrimination, KM and cohort tables. No plots."""
from pathlib import Path
import json, pickle, hashlib, importlib.metadata, os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from lifelines import KaplanMeierFitter
from lifelines.utils import median_survival_times
from survival_metrics import bootstrap_metrics, cindex, roc_auc

ROOT=Path(__file__).resolve().parents[1]
STUDY=Path(os.environ.get('HCC_OS_LEGACY_PROJECT_ROOT',ROOT.parents[1])).expanduser().resolve()
OUT=ROOT/'data/processed'
RAW=ROOT/'data/raw'
B=1000
SEED=20260905
TIMES=(365.,1095.,1825.)


def load_full(co):
    x=pd.read_csv(RAW/f'{co}_full_expression_logscale.csv.gz',index_col=0).T
    assert not x.index.duplicated().any() and not x.columns.duplicated().any()
    return x


def standardized_risk(model,x):
    prep=model['preprocessing'];z=x.loc[:,model['genes']].to_numpy(float)
    z[~np.isfinite(z)]=np.nan;z=np.clip(z,prep['lo'],prep['hi'])
    z=np.where(np.isnan(z),prep['mean'],z)
    keep=prep['valid_gene'];z=(z[:,keep]-prep['mean'][keep])/prep['sd'][keep]
    return model['model'].predict(z[:,prep['selected']])


def km_tables(d,co,predictor,groups):
    curves=[];risk=[]
    for group in ('Low','High'):
        a=d.loc[groups==group]
        if len(a)==0:continue
        km=KaplanMeierFitter().fit(a.time,event_observed=a.event,label='survival')
        tab=km.survival_function_.join(km.confidence_interval_).reset_index()
        tab.columns=['time_days','survival','ci_lower','ci_upper']
        tab['n_censor']=km.event_table.censored.to_numpy();tab['n_event']=km.event_table.observed.to_numpy();tab['n_risk']=km.event_table.at_risk.to_numpy()
        tab['cohort']=co;tab['predictor']=predictor;tab['group']=group;curves.append(tab)
        for t in np.arange(0,2191,365):
            risk.append(dict(cohort=co,predictor=predictor,group=group,time_days=t,n_risk=int((a.time>=t).sum()),n_total=len(a)))
    return pd.concat(curves),pd.DataFrame(risk)


def cohort_row(d,co):
    km=KaplanMeierFitter().fit(d.time, event_observed=1-d.event)
    medci=median_survival_times(km.confidence_interval_).iloc[0].to_numpy()
    return dict(cohort=co,n=len(d),events=int(d.event.sum()),censored=int((1-d.event).sum()),
      observed_time_min=d.time.min(),observed_time_max=d.time.max(),observed_time_median=d.time.median(),
      reverse_KM_followup_median_days=km.median_survival_time_,reverse_KM_followup_ci_lower_days=medci[0],reverse_KM_followup_ci_upper_days=medci[1],
      n_risk_1y=int((d.time>365).sum()),n_risk_3y=int((d.time>1095).sum()),n_risk_5y=int((d.time>1825).sum()))


def main():
    tc=load_full('TCGA');ic=load_full('ICGC')
    definitions=json.loads((ROOT/'references/signature_definitions.json').read_text())
    ros=set((ROOT/'references/HALLMARK_ROS_msigdb_v7.0_genes.txt').read_text().split())
    mp=pd.read_csv(ROOT.parent/'single_cell/data/processed/discovery_NMF_genes.csv')
    mp4=set(mp.loc[mp.program=='MP4','gene'])
    universe=sorted(set(tc.columns)&set(ic.columns))
    gene_sets={'ROS_mean_rank':ros,'MP4_mean_rank':mp4,'MP4_without_ROS_mean_rank':mp4-ros}
    mappings=[];corr=[];scores={};km=[];risk=[];cohort=[];cutoffs=[];availability=[]
    for label,gs in gene_sets.items():
        for g in sorted(gs):mappings.append(dict(gene_set=label,gene=g,TCGA_available=g in tc.columns,ICGC_available=g in ic.columns,used_in_common_universe=g in universe))
    pd.DataFrame(mappings).to_csv(OUT/'bulk_program_gene_coverage.csv',index=False)
    (OUT/'bulk_program_definition.json').write_text(json.dumps(dict(method='mean within-sample fractional expression rank (average ties) in a fixed shared gene universe',
        shared_universe_size=len(universe),ROS_gene_count=len(ros),MP4_gene_count=len(mp4),MP4_ROS_overlap=sorted(mp4&ros)),indent=2))
    (OUT/'bulk_shared_gene_universe.txt').write_text('\n'.join(universe)+'\n')
    with open(OUT/'new_tcga_only_final_model.pkl','rb') as f:new_model=pickle.load(f)
    for co,x in [('TCGA-LIHC',tc),('ICGC-LIRI',ic)]:
        d=pd.read_csv(RAW/f'{co}_frozen_risks.csv').rename(columns={'OS.time':'time','OS':'event','RS':'OSARS_frozen'}).set_index('ID')
        x=x.loc[d.index]
        assert np.all(np.isfinite(x)) and np.all(x>=0)
        signature_x=x.copy()
        if co=='ICGC-LIRI':
            recovered=pd.read_csv(RAW/'ICGC_unfiltered_published_signature_log_expression.csv',index_col=0)
            for gene in recovered.columns:signature_x[gene]=recovered.loc[d.index,gene]
        for label in ['Hong2022_8gene','Ma2024_TR_OSRG_3gene']:
            coef=pd.Series(definitions[label]['coefficients'])
            missing=sorted(set(coef.index)-set(signature_x.columns))
            availability.append(dict(cohort=co,signature=label,required_genes=len(coef),available_genes=len(coef)-len(missing),missing_genes=';'.join(missing)))
            if missing:continue
            raw=signature_x[coef.index]
            z=(raw-raw.mean())/raw.std(ddof=1)
            d[label+'_log']=raw@coef
            d[label+'_zscore']=z@coef
            # Frozen TCGA means/SD for transport sensitivity.
            d[label+'_TCGA_zscore']=((raw-tc[coef.index].mean())/tc[coef.index].std(ddof=1))@coef
        ranks=x[universe].rank(axis=1,method='average',pct=True)
        for label,gs in gene_sets.items():
            used=sorted(gs&set(universe));d[label]=ranks[used].mean(axis=1)
            r,p=spearmanr(d.OSARS_frozen,d[label]);
            # Patient/sample bootstrap CI for Spearman, fixed scores.
            rng=np.random.default_rng(SEED);v=[]
            for _ in range(B):
                ix=rng.integers(len(d),size=len(d));v.append(spearmanr(d.OSARS_frozen.iloc[ix],d[label].iloc[ix]).statistic)
            lo,hi=np.quantile(v,[.025,.975]);corr.append(dict(cohort=co,program=label,n=len(d),spearman_rho=r,p=p,ci_lower=lo,ci_upper=hi,genes_used=len(used)))
        d['new_TCGA_only_model']=standardized_risk(new_model,x)
        d['analysis_role']='training apparent' if co=='TCGA-LIHC' else 'historically consulted external cohort'
        scores[co]=d
        for label in ['OSARS_frozen','Hong2022_8gene_log','Ma2024_TR_OSRG_3gene_zscore']:
            if label not in d or d[label].isna().any():continue
            cutoff=d[label].median();g=np.where(d[label]>cutoff,'High','Low')
            d[label+'_group']=g
            cutoffs.append(dict(cohort=co,predictor=label,rule='score > cohort median',cutoff=cutoff,n_low=int((g=='Low').sum()),n_high=int((g=='High').sum())))
            k,r=km_tables(d,co,label,g);km.append(k);risk.append(r)
        # Transfer of original training cutoff: sensitivity, not optimized in validation.
        trcut=float(scores['TCGA-LIHC'].OSARS_frozen.median())
        g=np.where(d.OSARS_frozen>trcut,'High','Low');d['OSARS_training_cutoff_group']=g
        cutoffs.append(dict(cohort=co,predictor='OSARS_frozen_training_cutoff',rule='score > frozen TCGA median',cutoff=trcut,n_low=int((g=='Low').sum()),n_high=int((g=='High').sum())))
        k,r=km_tables(d,co,'OSARS_frozen_training_cutoff',g);km.append(k);risk.append(r)
        d.to_csv(OUT/f'{co}_all_scores.csv')
        cohort.append(cohort_row(d,co))
    # Unified benchmarking requires complete observed expression for BOTH formulas.
    icgc_benchmark=scores['ICGC-LIRI'].dropna(subset=['Hong2022_8gene_log','Ma2024_TR_OSRG_3gene_zscore']).copy()
    icgc_benchmark.to_csv(OUT/'ICGC-LIRI_benchmark_completecases_all_scores.csv')
    pd.DataFrame({'ID':scores['ICGC-LIRI'].index,'complete_for_both_signatures':scores['ICGC-LIRI'].index.isin(icgc_benchmark.index)}).to_csv(OUT/'ICGC_benchmark_inclusion.csv',index=False)
    scores['ICGC-LIRI_benchmark_completecases']=icgc_benchmark
    cohort.append(cohort_row(icgc_benchmark,'ICGC-LIRI_benchmark_completecases'))
    # Exact expression duplicates are detectable; donor identity is unavailable.
    hashes=pd.util.hash_pandas_object(ic,index=False)
    dup=pd.DataFrame({'ID':ic.index,'expression_profile_hash':hashes.to_numpy()})
    dup['identical_expression_profile_n']=dup.expression_profile_hash.map(dup.expression_profile_hash.value_counts())
    dup['retain_profile_sensitivity']=~dup.expression_profile_hash.duplicated()
    dup.to_csv(OUT/'ICGC_identical_expression_profile_audit.csv',index=False)
    keep=dup.loc[dup.retain_profile_sensitivity,'ID']
    scores['ICGC-LIRI_unique_profiles']=scores['ICGC-LIRI'].loc[keep].copy()
    scores['ICGC-LIRI_unique_profiles'].to_csv(OUT/'ICGC-LIRI_unique_profiles_all_scores.csv')
    cohort.append(cohort_row(scores['ICGC-LIRI_unique_profiles'],'ICGC-LIRI_unique_profiles'))
    # CPTAC surrogate: frozen correlation-weighted protein score, NOT GBM input.
    c=pd.read_csv(STUDY/'CPTAC_proteomics_OSARS/CPTAC_OSARS_protein_surrogate_scores_OS_RFS.csv')
    for endpoint,tcol,ecol in [('OS','ostime','osevent'),('RFS','rsftime','rsfevent')]:
        co='CPTAC_'+endpoint
        d=pd.DataFrame({'ID':c.Sample,'time':c[tcol]*365/12,'event':c[ecol].astype(int),'OSARS_protein_surrogate':c.OSARS_protein}).set_index('ID')
        scores[co]=d;g=np.where(d.OSARS_protein_surrogate>d.OSARS_protein_surrogate.median(),'High','Low')
        d['OSARS_protein_surrogate_group']=g;d.to_csv(OUT/f'{co}_all_scores.csv')
        cutoffs.append(dict(cohort=co,predictor='OSARS_protein_surrogate',rule='score > cohort median',cutoff=d.OSARS_protein_surrogate.median(),n_low=int((g=='Low').sum()),n_high=int((g=='High').sum())))
        k,r=km_tables(d,co,'OSARS_protein_surrogate',g);km.append(k);risk.append(r);cohort.append(cohort_row(d,co))
    # Each cohort uses exactly the same eligible rows for all compared markers.
    metrics=[];draws=[];coords=[];diffs=[]
    for co,d in scores.items():
        cols=['OSARS_protein_surrogate'] if co.startswith('CPTAC') else ['OSARS_frozen','Hong2022_8gene_log','Ma2024_TR_OSRG_3gene_zscore','ROS_mean_rank','MP4_mean_rank','MP4_without_ROS_mean_rank','new_TCGA_only_model','Hong2022_8gene_zscore','Ma2024_TR_OSRG_3gene_log']
        cols=[x for x in cols if x in d and d[x].notna().all()]
        print('Bootstrap',co,len(d),len(cols),flush=True)
        a,b,cube=bootstrap_metrics(d.time,d.event,{k:d[k] for k in cols},TIMES,B,SEED)
        a['cohort']=co;b['cohort']=co;metrics.append(a);draws.append(b)
        for label in cols:
            for horizon in TIMES:
                _,r=roc_auc(d.time,d.event,d[label],horizon,True)
                if r is not None:r['cohort']=co;r['predictor']=label;r['time_days']=horizon;coords.append(r)
        if 'OSARS_frozen' in cols:
            for j,label in enumerate(cols[1:],1):
                for k,(metric,t) in enumerate([('C-index',np.nan)]+[('AUC',t) for t in TIMES]):
                    delta=cube[:,0,k]-cube[:,j,k];good=delta[np.isfinite(delta)]
                    lo,hi=np.quantile(good,[.025,.975]) if len(good) else [np.nan,np.nan]
                    aa=a[(a.metric==metric)&(a.predictor=='OSARS_frozen')]
                    bb=a[(a.metric==metric)&(a.predictor==label)]
                    if metric=='AUC':aa=aa[aa.time_days==t];bb=bb[bb.time_days==t]
                    diffs.append(dict(cohort=co,reference='OSARS_frozen',comparator=label,metric=metric,time_days=t,
                      delta=float(aa.estimate.iloc[0]-bb.estimate.iloc[0]),ci_lower=lo,ci_upper=hi,bootstrap_valid=len(good),
                      interpretation='descriptive paired contrast; model-selection/development overlap not removed'))
    pd.concat(metrics).to_csv(OUT/'discrimination_metrics_with_CI.csv',index=False)
    pd.concat(draws).to_csv(OUT/'discrimination_bootstrap_draws.csv.gz',index=False)
    pd.concat(coords).to_csv(OUT/'ROC_coordinates.csv',index=False)
    pd.concat(km).to_csv(OUT/'KM_coordinates.csv',index=False)
    pd.concat(risk).to_csv(OUT/'KM_risk_tables.csv',index=False)
    pd.DataFrame(cohort).to_csv(OUT/'cohort_survival_summary.csv',index=False)
    pd.DataFrame(cutoffs).to_csv(OUT/'risk_cutoff_manifest.csv',index=False)
    pd.DataFrame(diffs).to_csv(OUT/'paired_discrimination_differences.csv',index=False)
    pd.DataFrame(corr).to_csv(OUT/'OSARS_program_correlations.csv',index=False)
    pd.DataFrame(availability).to_csv(OUT/'published_signature_gene_coverage.csv',index=False)
    # Closed-loop comparison with the previously generated timeROC results.
    old=pd.read_csv(STUDY/'FIG3_ROC_Cindex/fig3_timeROC_AUC_summary.csv')
    a=pd.concat(metrics);a=a[(a.predictor=='OSARS_frozen')&(a.metric=='AUC')]
    check=old.merge(a,left_on=['dataset','time_days'],right_on=['cohort','time_days'])
    check['absolute_error']=abs(check.AUC-check.estimate)
    assert check.absolute_error.max()<1e-12
    check.to_csv(OUT/'timeROC_point_estimate_reproduction_check.csv',index=False)
    packages=['numpy','pandas','scipy','lifelines','scikit-learn','scikit-survival']
    (OUT/'python_versions.json').write_text(json.dumps({p:importlib.metadata.version(p) for p in packages},indent=2))


if __name__=='__main__':main()
