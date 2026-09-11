"""Bounded MP4 confounding sensitivities and verified PDC tumor-only survival.
No plot generation; raw and previously processed source files are read only.
"""
from pathlib import Path
from importlib import import_module
import numpy as np
import pandas as pd
from scipy import stats
from lifelines.statistics import logrank_test
from survival_metrics import bootstrap_metrics, roc_auc

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/processed'
B=1000;SEED=20260905
scores=import_module('05_prepare_scores_and_metrics')

def partial_rho(d,program,continuous,stage):
    x=stats.rankdata(d.OSARS_frozen);y=stats.rankdata(d[program])
    z=[np.ones(len(d))]+[stats.rankdata(d[v]) for v in continuous]
    if stage:
        z.extend(pd.get_dummies(d.stage_num.astype(int),drop_first=True).to_numpy(float).T)
    z=np.asarray(z).T
    rx=x-z@np.linalg.lstsq(z,x,rcond=None)[0]
    ry=y-z@np.linalg.lstsq(z,y,rcond=None)[0]
    r=np.corrcoef(rx,ry)[0,1];k=np.linalg.matrix_rank(z)-1
    p=2*stats.t.sf(abs(r)*np.sqrt((len(d)-k-2)/(1-r*r)),len(d)-k-2)
    return r,p,k,rx,ry

def partial_analysis():
    purity=pd.read_csv(ROOT.parent/'immune_therapy/data/processed/immune_patient_values.csv')
    purity=purity.loc[purity.Feature=='Tumor purity',['Sample','Value']].rename(columns={'Sample':'ID','Value':'TumorPurity'})
    rows=[];residuals=[]
    for co in ['TCGA-LIHC','ICGC-LIRI']:
        d=pd.read_csv(OUT/f'{co}_clinical_analysis_rows.csv')
        if co=='TCGA-LIHC':d=d.merge(purity,on='ID',how='left',validate='one_to_one')
        specs=[('unadjusted',[],False),('ROS',['ROS_mean_rank'],False),
               ('ROS + stage',['ROS_mean_rank'],True)]
        if co=='TCGA-LIHC':specs.append(('ROS + stage + ESTIMATE purity',['ROS_mean_rank','TumorPurity'],True))
        for program in ['MP4_mean_rank','MP4_without_ROS_mean_rank']:
            for label,cont,stage in specs:
                needed=['OSARS_frozen',program]+cont+(['stage_num'] if stage else [])
                a=d.dropna(subset=needed).copy()
                r,p,k,rx,ry=partial_rho(a,program,cont,stage)
                rng=np.random.default_rng(SEED);boot=[]
                for _ in range(B):
                    ix=rng.integers(len(a),size=len(a));boot.append(partial_rho(a.iloc[ix],program,cont,stage)[0])
                lo,hi=np.nanquantile(boot,[.025,.975])
                rows.append(dict(cohort=co,program=program,adjustment=label,n=len(a),covariate_df=k,partial_spearman_rho=r,p_approximate=p,ci_lower=lo,ci_upper=hi,bootstrap_B=B,
                    method='Pearson correlation of rank residuals; continuous covariates ranked, stage categorical; two-sided t approximation; sample bootstrap reranks and refits'))
                b=pd.DataFrame(dict(ID=a.ID,OSARS_rank_residual=rx,program_rank_residual=ry))
                b['cohort']=co;b['program']=program;b['adjustment']=label;residuals.append(b)
    tab=pd.DataFrame(rows);tab['q_BH_all_14_tests']=stats.false_discovery_control(tab.p_approximate)
    tab.to_csv(OUT/'MP4_partial_rank_correlations.csv',index=False)
    pd.concat(residuals).to_csv(OUT/'MP4_partial_rank_residual_coordinates.csv',index=False)

def pdc_analysis():
    # Metadata discrepancy is independently verified by the protein agent.
    audit=pd.read_csv(ROOT.parent/'protein_support/data/processed/pdc_source_tissue_discrepancies.csv')
    exclusions=audit.loc[(audit.sample_id.str.startswith('T'))&(audit.sample_type!='Primary Tumor'),'sample_id'].tolist()
    assert exclusions==['T724']
    mets=[];curves=[];risks=[];rocs=[];summaries=[];cutoffs=[];tests=[]
    for ep in ['OS','RFS']:
        d=pd.read_csv(OUT/f'CPTAC_{ep}_clinical_analysis_rows.csv')
        d=d.loc[~d.ID.isin(exclusions)].copy();assert len(d)==158
        co=f'PDC000198_{ep}_verified158'
        d['OSARS_SD']=(d.OSARS_protein_surrogate-d.OSARS_protein_surrogate.mean())/d.OSARS_protein_surrogate.std(ddof=1)
        cutoff=d.OSARS_protein_surrogate.median();d['High']=(d.OSARS_protein_surrogate>cutoff).astype(int)
        d['OSARS_protein_surrogate_group']=np.where(d.High,'High','Low')
        d.to_csv(OUT/f'{co}_all_scores.csv',index=False)
        s,_,_=bootstrap_metrics(d.time,d.event,{'OSARS_protein_surrogate':d.OSARS_protein_surrogate},(365.,1095.),B,SEED)
        s['cohort']=co;mets.append(s)
        for t in [365.,1095.]:
            _,r=roc_auc(d.time,d.event,d.OSARS_protein_surrogate,t,True)
            if r is not None:r['time_days']=t;r['cohort']=co;r['predictor']='OSARS_protein_surrogate';rocs.append(r)
        k,r=scores.km_tables(d,co,'OSARS_protein_surrogate',d.OSARS_protein_surrogate_group)
        curves.append(k);risks.append(r);summaries.append(scores.cohort_row(d,co))
        lo=d.loc[d.High==0];hi=d.loc[d.High==1];test=logrank_test(lo.time,hi.time,lo.event,hi.event)
        tests.append(dict(cohort=co,grouping='High',n=len(d),events=int(d.event.sum()),logrank_chisq=test.test_statistic,df=1,p=test.p_value))
        cutoffs.append(dict(cohort=co,predictor='OSARS_protein_surrogate',rule='score > cohort median after verified-tumor filter',cutoff=cutoff,n_low=int((d.High==0).sum()),n_high=int((d.High==1).sum())))
    for name,frames in [('PDC_verified158_discrimination_CI',mets),('PDC_verified158_KM_coordinates',curves),('PDC_verified158_KM_risk_tables',risks),('PDC_verified158_ROC_coordinates',rocs)]:
        pd.concat(frames).to_csv(OUT/f'{name}.csv',index=False)
    for name,values in [('PDC_verified158_cohort_survival_summary',summaries),('PDC_verified158_cutoff_manifest',cutoffs),('PDC_verified158_KM_logrank_tests',tests)]:
        pd.DataFrame(values).to_csv(OUT/f'{name}.csv',index=False)
    pd.DataFrame({'excluded_archived_tumor_label':exclusions,'reason':'PDC biospecimen metadata identifies solid tissue normal'}).to_csv(OUT/'PDC_verified158_exclusions.csv',index=False)

if __name__=='__main__':
    partial_analysis();pdc_analysis()
