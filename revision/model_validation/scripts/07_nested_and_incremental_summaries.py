"""Conditional bootstrap intervals for OOF discrimination and clinical increments."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from survival_metrics import bootstrap_metrics,cindex,roc_auc
from importlib import import_module
scores=import_module('05_prepare_scores_and_metrics')
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/processed'
B=1000;SEED=20260905


def main():
    d=pd.read_csv(OUT/'nested_tcga_oof_predictions.csv')
    folds=pd.read_csv(OUT/'nested_outer_sample_splits.csv');folds=folds[folds.role=='test']
    d=d.merge(folds[['ID','outer_fold']],on='ID',validate='one_to_one')
    frozen=pd.read_csv(OUT/'TCGA-LIHC_all_scores.csv')[['ID','OSARS_frozen']]
    d=d.merge(frozen,on='ID',validate='one_to_one')
    # Primary aggregate: equal mean within-fold C, since different fold models
    # have different link scales. Resample held-out subjects within each fold.
    means=[];rows=[];old_means=[]
    for fold,a in d.groupby('outer_fold'):
        a=a.reset_index(drop=True)
        s,b,cube=bootstrap_metrics(a.time,a.event,{'OOF':a.OOF_link_score,'Frozen_apparent':a.OSARS_frozen},(),B,SEED+int(fold))
        s['outer_fold']=fold;rows.append(s);means.append(cube[:,0,0]);old_means.append(cube[:,1,0])
    samples=np.mean(means,axis=0);old_samples=np.mean(old_means,axis=0)
    primary=float(pd.concat(rows).query("predictor=='OOF'").estimate.mean())
    apparent=float(pd.concat(rows).query("predictor=='Frozen_apparent'").estimate.mean())
    summary=pd.DataFrame([
      dict(predictor='Nested CV: mean outer-fold C-index',estimate=primary,ci_lower=np.quantile(samples,.025),ci_upper=np.quantile(samples,.975),n=len(d),events=int(d.event.sum()),interpretation='conditional on fixed outer fits and independently derived candidate gene set; excludes resampling model-training uncertainty'),
      dict(predictor='Frozen legacy: same fold-wise apparent C-index',estimate=apparent,ci_lower=np.quantile(old_samples,.025),ci_upper=np.quantile(old_samples,.975),n=len(d),events=int(d.event.sum()),interpretation='training apparent; comparison changes training strategy, not an unbiased optimism correction of the frozen model'),
      dict(predictor='Apparent minus nested training-strategy contrast',estimate=apparent-primary,ci_lower=np.quantile(old_samples-samples,.025),ci_upper=np.quantile(old_samples-samples,.975),n=len(d),events=int(d.event.sum()),interpretation='descriptive contrast of different fitting strategies; not estimator of pure optimism for one fixed model')])
    summary.to_csv(OUT/'nested_Cindex_summary_CI.csv',index=False)
    pd.concat(rows).to_csv(OUT/'nested_outer_fold_Cindex_CI.csv',index=False)
    # Pooled ROC/KM uses training-percentile normalized scores to align fold scales.
    a,b,_=bootstrap_metrics(d.time,d.event,{'OOF_training_percentile':d.OOF_training_percentile},(365.,1095.,1825.),B,SEED)
    a['interpretation']='pooled descriptive OOF score normalized by each training-fold empirical CDF; primary C-index is mean within-fold'
    a.to_csv(OUT/'nested_pooled_percentile_discrimination_CI.csv',index=False)
    co=[]
    for h in (365.,1095.,1825.):
        _,r=roc_auc(d.time,d.event,d.OOF_training_percentile,h,True)
        if r is not None:r['time_days']=h;co.append(r)
    pd.concat(co).to_csv(OUT/'nested_OOF_ROC_coordinates.csv',index=False)
    k,r=scores.km_tables(d,'TCGA_nested_OOF','OOF_training_cutoff',d.OOF_group)
    k.to_csv(OUT/'nested_OOF_KM_coordinates.csv',index=False);r.to_csv(OUT/'nested_OOF_KM_risk_tables.csv',index=False)
    # Stage increments use identical complete cases per cohort; all coefficients
    # came from TCGA only (R script 06). These still concern the frozen legacy score.
    d=pd.read_csv(OUT/'clinical_incremental_predictions.csv');met=[];diff=[];cur=[]
    for co,a in d.groupby('cohort'):
        wide=a.pivot(index=['ID','time','event'],columns='predictor',values='score').reset_index()
        cols=['Stage','Stage_plus_ROS','Stage_plus_OSARS']
        s,b,cube=bootstrap_metrics(wide.time,wide.event,{p:wide[p] for p in cols},(365.,1095.,1825.),B,SEED)
        s['cohort']=co;met.append(s)
        for j in [1,2]:
            for k,(metric,h) in enumerate([('C-index',np.nan),('AUC',365.),('AUC',1095.),('AUC',1825.)]):
                v=cube[:,j,k]-cube[:,0,k];v=v[np.isfinite(v)];lo,hi=np.quantile(v,[.025,.975]) if len(v) else [np.nan,np.nan]
                sa=s[(s.predictor==cols[j])&(s.metric==metric)];sb=s[(s.predictor=='Stage')&(s.metric==metric)]
                if metric=='AUC':sa=sa[sa.time_days==h];sb=sb[sb.time_days==h]
                diff.append(dict(cohort=co,predictor=cols[j],reference='Stage',metric=metric,time_days=h,delta=float(sa.estimate.iloc[0]-sb.estimate.iloc[0]),ci_lower=lo,ci_upper=hi,bootstrap_valid=len(v)))
        for p in cols:
            for h in (365.,1095.,1825.):
                _,r=roc_auc(wide.time,wide.event,wide[p],h,True)
                if r is not None:r['cohort']=co;r['predictor']=p;r['time_days']=h;cur.append(r)
    pd.concat(met).to_csv(OUT/'clinical_incremental_discrimination_CI.csv',index=False)
    pd.DataFrame(diff).to_csv(OUT/'clinical_incremental_paired_differences.csv',index=False)
    pd.concat(cur).to_csv(OUT/'clinical_incremental_ROC_coordinates.csv',index=False)
    cfg=json.loads((OUT/'nested_cv_protocol.json').read_text())
    cfg['candidate_origin']='657 positive OS-high scRNA markers intersected with genes available in TCGA, ICGC and GSE14520 exactly reproduces ordered 448 vector; no Scissor/outcome filter identified in this step'
    (OUT/'nested_cv_protocol.json').write_text(json.dumps(cfg,indent=2))


if __name__=='__main__':main()
