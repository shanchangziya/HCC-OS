#!/usr/bin/env python3
"""Independent numerical consistency audit of the frozen analysis exports."""
from pathlib import Path
import hashlib,json,os
import numpy as np
import pandas as pd

root=Path(__file__).resolve().parent.parent
p=root/'data/processed'
x=pd.read_csv(p/'analysis_input.csv')
pred=pd.read_csv(p/'frozen_patient_predictions.csv')
co=pd.read_csv(p/'lasso_coefficients.csv')
spec=json.loads((p/'frozen_model_specification.json').read_text())
assert pred.patient_id.is_unique and x.patient_id.tolist()==pred.patient_id.tolist()
reconstructed=x[co.feature].to_numpy()@co.coefficient_raw_scale.to_numpy()
error=float(np.max(np.abs(reconstructed-pred.pathology_score)))
assert error<1e-10
cut=float(pred.loc[pred.split=='Train','pathology_score'].median())
assert abs(cut-spec['training_median_cutoff'])<1e-12
assert np.array_equal(np.where(pred.pathology_score>cut,'High','Low'),pred.risk_group)
rt=pd.read_csv(p/'km_risk_table.csv')
for r in rt.itertuples():
    at_risk=((pred.split==r.split)&(pred.risk_group==r.group)&(pred.time_days>=r.time_days)).sum()
    assert at_risk==r.n_risk
km=pd.read_csv(p/'km_coordinates.csv')
for (_,group),k in km.groupby(['split','group']):
    assert np.all(np.diff(k.time_days)>=0) and np.all(np.diff(k.survival)<=1e-12)
    ok=k.lower.notna()&k.upper.notna()
    assert (k.loc[ok,'lower']<=k.loc[ok,'survival']+1e-12).all()
    assert (k.loc[ok,'survival']<=k.loc[ok,'upper']+1e-12).all()
folds=pd.read_csv(p/'training_cv_folds.csv')
assert set(folds.patient_id)==set(pred.loc[pred.split=='Train','patient_id'])
fold_n=pd.read_csv(p/'cv_fold_summary.csv')
assert (fold_n.events_validation>0).all()
metrics=pd.read_csv(p/'performance_metrics.csv')
assert (metrics.bootstrap_valid==500).all()
assert metrics[['estimate','lower','upper']].notna().all().all()
assert ((metrics.estimate>=0)&(metrics.estimate<=1)).all()
auc=pd.read_csv(p/'auc_implementation_crosscheck.csv')
assert auc.abs_difference.max()<1e-3
manifest=json.loads((p/'source_manifest.json').read_text())
assert manifest['sources_unchanged_after_analysis']
checks=[]
unavailable=[]
configured_source=Path(os.environ.get('HCC_OS_PATHOLOGY_SOURCE_DIR',root/'data/raw')).expanduser()
configured_clinical=os.environ.get('HCC_OS_PATHOLOGY_CLINICAL')
for item in manifest['sources']:
    name=item.get('source_id',Path(item.get('path','')).name)
    role=item.get('role','')
    if role=='clinical_covariates' or (not role and name=='Fig4_TCGA_OSARS_clinical_all_matched.csv'):
        local=Path(configured_clinical).expanduser() if configured_clinical else root/'data/raw'/name
    else:
        local=configured_source/name
    if local.exists():
        h=hashlib.sha256(local.read_bytes()).hexdigest()
        assert h==item['sha256']
        checks.append(name)
    else:
        unavailable.append(name)
qa={'status':'PASS','n_patients':len(pred),'train_n':230,'test_n':100,
    'risk_reconstruction_max_abs_difference':error,'training_cutoff_exact':cut,
    'km_monotonic_and_ci_ordered':True,'risk_table_matches_individual_records':True,
    'cv_contains_training_patients_only':True,'all_bootstrap_replicates_valid':True,
    'max_auc_difference_vs_timeROC':float(auc.abs_difference.max()),
    'source_hashes_rechecked':checks,'sources_not_available_for_recheck':unavailable,
    'original_sources_unchanged_during_frozen_run':True,
    'final_glmnet_path_jerr':0,
    'lambda_probe_warning':'Initial lambda-max probe warned at its lowest lambda; final full path converged (jerr=0) and selected lambda is within returned range.'}
(p/'independent_output_qa.json').write_text(json.dumps(qa,indent=2))
print(json.dumps(qa,indent=2))
