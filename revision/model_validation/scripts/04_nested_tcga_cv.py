"""TCGA-only nested CV conditional on the archived scRNA candidate gene set.

This is a NEW, bounded training strategy, not validation of the frozen GBM and
not a re-run of all 117 legacy configurations. No ICGC outcomes are read here.
No plots are generated. All splits, preprocessing, selected genes, tuning
results and held-out scores are exported.
"""
from pathlib import Path
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json, pickle, time
import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.model_selection import StratifiedKFold
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sksurv.ensemble import GradientBoostingSurvivalAnalysis
from sksurv.util import Surv
from survival_metrics import cindex, bootstrap_metrics

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/processed'
SEED=20260905


def univariate_cox_pvalues(x, t, event):
    """Simultaneous independent one-variable Cox fits with Efron ties.

    x is already scaled using the current training sample only. Uses a Newton
    iteration with bounded steps; marks failed fits, never labels them selected.
    """
    order=np.argsort(t,kind='stable'); x=x[order];t=t[order];event=event[order]
    unique=np.unique(t[event]); risk_index=np.searchsorted(t,unique)
    event_indices=[np.flatnonzero((t==v)&event) for v in unique]
    beta=np.zeros(x.shape[1]); converged=np.zeros(x.shape[1],bool)
    for iteration in range(60):
        eta=np.clip(x*beta,-50,50);w=np.exp(eta)
        s0=np.cumsum(w[::-1],axis=0)[::-1][risk_index]
        s1=np.cumsum((w*x)[::-1],axis=0)[::-1][risk_index]
        s2=np.cumsum((w*x*x)[::-1],axis=0)[::-1][risk_index]
        u=np.zeros(x.shape[1]);info=np.zeros(x.shape[1])
        for k,ix in enumerate(event_indices):
            d=len(ix);d0=w[ix].sum(axis=0);d1=(w[ix]*x[ix]).sum(axis=0);d2=(w[ix]*x[ix]**2).sum(axis=0)
            u+=x[ix].sum(axis=0)
            for f in np.arange(d)/d:
                z0=s0[k]-f*d0;z1=s1[k]-f*d1;z2=s2[k]-f*d2
                u-=z1/z0;info+=z2/z0-(z1/z0)**2
        step=np.divide(u,info,out=np.zeros_like(u),where=info>1e-10)
        converged=np.abs(step)<1e-7
        if converged.all():break
        beta+=np.clip(step,-2,2)
    z=beta*np.sqrt(np.maximum(info,0))
    p=2*norm.sf(np.abs(z));p[~converged]=np.nan
    return beta,p,converged


def prepare(train, valid, t, e):
    tr=np.asarray(train,float).copy();va=np.asarray(valid,float).copy()
    tr[~np.isfinite(tr)]=np.nan;va[~np.isfinite(va)]=np.nan
    lo,hi=np.nanquantile(tr,[.01,.99],axis=0)
    tr=np.clip(tr,lo,hi);va=np.clip(va,lo,hi)
    mean=np.nanmean(tr,axis=0)
    tr=np.where(np.isnan(tr),mean,tr);va=np.where(np.isnan(va),mean,va)
    sd=tr.std(axis=0,ddof=1);valid_gene=np.isfinite(sd)&(sd>1e-6)
    tr=(tr[:,valid_gene]-mean[valid_gene])/sd[valid_gene]
    va=(va[:,valid_gene]-mean[valid_gene])/sd[valid_gene]
    beta,p,converged=univariate_cox_pvalues(tr,t,e)
    selected=np.isfinite(p)&(p<.01)
    info={'lo':lo,'hi':hi,'mean':mean,'sd':sd,'valid_gene':valid_gene,'selected':selected,'beta':beta,'p':p,'converged':converged}
    return tr[:,selected],va[:,selected],info


GRID=[
 {'name':'Ridge_Cox_penalty1','family':'ridge','penalty':1.},
 {'name':'Ridge_Cox_penalty10','family':'ridge','penalty':10.},
 {'name':'ElasticNet_Cox_penalty0.01','family':'elastic','penalty':.01},
 {'name':'ElasticNet_Cox_penalty0.1','family':'elastic','penalty':.1},
 *[{'name':f'Cox_GBM_depth{d}_trees{n}','family':'gbm','depth':d,'n':n} for d in (1,3) for n in (100,400)]
]


def fit_model(spec,x,y,seed):
    if x.shape[1]==0:return None
    if spec['family']=='ridge':
        model=CoxPHSurvivalAnalysis(alpha=spec['penalty'],ties='efron',n_iter=200)
    elif spec['family']=='elastic':
        model=CoxnetSurvivalAnalysis(l1_ratio=.5,alphas=[spec['penalty']],max_iter=100000,normalize=False)
    else:
        model=GradientBoostingSurvivalAnalysis(loss='coxph',learning_rate=.03,n_estimators=spec['n'],
          max_depth=spec['depth'],min_samples_leaf=10,subsample=.5,random_state=seed)
    model.fit(x,y);return model


def prediction(model,x):
    return np.zeros(len(x)) if model is None else model.predict(x)


def tune(x,t,e,outer_id):
    cv=StratifiedKFold(n_splits=3,shuffle=True,random_state=SEED+outer_id)
    scores=np.full((3,len(GRID)),np.nan);records=[];splits=[]
    for fold,(tr,va) in enumerate(cv.split(x,e),1):
        a,b,prep=prepare(x[tr],x[va],t[tr],e[tr]); y=Surv.from_arrays(e[tr],t[tr])
        for j,spec in enumerate(GRID):
            err=''
            try:
                model=fit_model(spec,a,y,SEED+outer_id*100+fold)
                scores[fold-1,j]=cindex(t[va],e[va],prediction(model,b))
            except Exception as exc:err=str(exc)
            records.append(dict(outer_fold=outer_id,inner_fold=fold,model=spec['name'],cindex=scores[fold-1,j],n_train=len(tr),n_validate=len(va),selected_genes=a.shape[1],error=err))
        splits.extend(dict(outer_fold=outer_id,inner_fold=fold,training_row=int(i),role='train') for i in tr)
        splits.extend(dict(outer_fold=outer_id,inner_fold=fold,training_row=int(i),role='validation') for i in va)
    means=np.mean(scores,axis=0) # any failure disqualifies the configuration
    chosen=int(np.nanargmax(means)) # ties resolve by predeclared GRID order
    return chosen,records,splits


def main():
    tic=time.time();OUT.mkdir(parents=True,exist_ok=True)
    xdf=pd.read_csv(ROOT/'data/raw/TCGA_unclipped_448_expression.csv').set_index('ID')
    dat=pd.read_csv(ROOT/'data/raw/TCGA-LIHC_frozen_risks.csv').set_index('ID')
    assert not xdf.index.duplicated().any() and set(xdf.index)==set(dat.index)
    x=xdf.loc[dat.index].to_numpy(float);genes=xdf.columns.to_numpy()
    t=dat['OS.time'].to_numpy(float);e=dat.OS.to_numpy(bool)
    assert np.all(t>0)
    cfg=dict(seed=SEED,outer_folds=5,inner_folds=3,selection='training-only univariate Efron Cox Wald p < 0.01',
      preprocessing='current training 1%/99% quantile winsorization; training mean imputation; training mean/SD standardization',
      criterion='mean of three inner-validation Harrell C-indices',grid=GRID,
      scope='NEW bounded training strategy; conditional on archived 448 candidate set; not re-run of legacy 117 model search',
      external_outcomes_used=False,candidate_origin='657 positive OS-high scRNA markers intersected with genes available in TCGA, ICGC and GSE14520 exactly reproduces ordered 448 vector; no Scissor/outcome filter identified in this step')
    (OUT/'nested_cv_protocol.json').write_text(json.dumps(cfg,indent=2))
    folds=StratifiedKFold(n_splits=5,shuffle=True,random_state=SEED)
    oof=np.full(len(x),np.nan);oof_percentile=np.full(len(x),np.nan);oof_group=np.full(len(x),'',object)
    outer_rows=[];inner_rows=[];inner_splits=[];sample_splits=[];selection=[];preprocess=[]
    for outer,(tr,te) in enumerate(folds.split(x,e),1):
        print('Outer fold',outer,'start',round(time.time()-tic,1),flush=True)
        j,rows,split=tune(x[tr],t[tr],e[tr],outer);inner_rows+=rows
        for r in split:r['ID']=str(dat.index[tr[r.pop('training_row')]])
        inner_splits+=split
        a,b,prep=prepare(x[tr],x[te],t[tr],e[tr]);spec=GRID[j]
        model=fit_model(spec,a,Surv.from_arrays(e[tr],t[tr]),SEED+outer*1000)
        pred=prediction(model,b);predtr=prediction(model,a)
        cutoff=float(np.median(predtr)); oof[te]=pred
        # Convert heterogeneous outer models to training-distribution percentiles
        # for pooled descriptive KM/ROC; primary discrimination is within-fold.
        oof_percentile[te]=np.searchsorted(np.sort(predtr),pred,side='right')/len(predtr)
        oof_group[te]=np.where(pred>cutoff,'High','Low')
        outer_rows.append(dict(outer_fold=outer,selected_model=spec['name'],n_train=len(tr),n_test=len(te),test_events=int(e[te].sum()),selected_genes=a.shape[1],cindex=cindex(t[te],e[te],pred),training_median_cutoff=cutoff))
        sample_splits.extend(dict(ID=dat.index[i],outer_fold=outer,role='train') for i in tr)
        sample_splits.extend(dict(ID=dat.index[i],outer_fold=outer,role='test') for i in te)
        for g,be,pv,conv,sel in zip(genes[prep['valid_gene']],prep['beta'],prep['p'],prep['converged'],prep['selected']):
            selection.append(dict(outer_fold=outer,gene=g,beta=be,p=pv,converged=conv,selected=sel))
        for g,lo,hi,mu,sd in zip(genes,prep['lo'],prep['hi'],prep['mean'],prep['sd']):
            preprocess.append(dict(outer_fold=outer,gene=g,q01=lo,q99=hi,mean=mu,sd=sd))
        with open(OUT/f'nested_outer{outer}_fit.pkl','wb') as f:pickle.dump(dict(model=model,preprocessing=prep,spec=spec,genes=genes),f)
        pd.DataFrame(outer_rows).to_csv(OUT/'nested_outer_fold_performance.csv',index=False)
    out=pd.DataFrame(dict(ID=dat.index,time=t,event=e.astype(int),OOF_link_score=oof,OOF_training_percentile=oof_percentile,OOF_group=oof_group))
    out.to_csv(OUT/'nested_tcga_oof_predictions.csv',index=False)
    pd.DataFrame(inner_rows).to_csv(OUT/'nested_inner_model_selection.csv',index=False)
    pd.DataFrame(inner_splits).to_csv(OUT/'nested_inner_sample_splits.csv',index=False)
    pd.DataFrame(sample_splits).to_csv(OUT/'nested_outer_sample_splits.csv',index=False)
    pd.DataFrame(selection).to_csv(OUT/'nested_outer_gene_selection.csv',index=False)
    pd.DataFrame(preprocess).to_csv(OUT/'nested_outer_preprocessing.csv',index=False)
    # Final strategy selection/training uses full TCGA only, AFTER OOF evaluation.
    j,rows,_=tune(x,t,e,100);a,_,prep=prepare(x,x,t,e)
    final=fit_model(GRID[j],a,Surv.from_arrays(e,t),SEED+100000)
    final_dict=dict(model=final,preprocessing=prep,spec=GRID[j],genes=genes,training_ids=dat.index.to_list(),training_predictions=prediction(final,a))
    with open(OUT/'new_tcga_only_final_model.pkl','wb') as f:pickle.dump(final_dict,f)
    pd.DataFrame(rows).to_csv(OUT/'new_final_inner_model_selection.csv',index=False)
    cfg['final_selected_model']=GRID[j];cfg['outer_mean_cindex']=float(pd.DataFrame(outer_rows).cindex.mean());cfg['elapsed_seconds']=time.time()-tic
    (OUT/'nested_cv_protocol.json').write_text(json.dumps(cfg,indent=2))
    print(json.dumps(cfg,indent=2),flush=True)


if __name__=='__main__':main()
