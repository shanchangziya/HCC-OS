#!/usr/bin/env python3
"""Verify pairing, exclusions, score arithmetic, and numerical report consistency."""
from pathlib import Path
import hashlib,json,math
import numpy as np
import pandas as pd

root=Path(__file__).resolve().parent.parent
p=root/'data/processed'
pairs=pd.read_csv(p/'nqo1_verified_pairs.csv')
all_samples=pd.read_csv(p/'nqo1_all_samples_with_verified_mapping.csv')
assert pairs.case_id.is_unique and len(pairs)==159
matched=all_samples[all_samples.case_id.isin(pairs.case_id)]
assert len(matched)==318 and matched.tissue_label_consistent.all()
assert (matched.groupby('case_id').tissue.nunique()==2).all()
assert matched.originally_measured.all()
assert np.allclose(pairs.tumor-pairs.adjacent,pairs.difference_tumor_minus_adjacent)
exclusions=pd.read_csv(p/'nqo1_unpaired_excluded_from_paired_test.csv')
assert len(exclusions)==12 and set(exclusions.sample_id).isdisjoint(set(matched.sample_id))
conflicts=pd.read_csv(p/'pdc_source_tissue_discrepancies.csv')
assert set(conflicts.sample_id)=={'T724','P723'}
surv=pd.read_csv(p/'nqo1_survival_input.csv')
verified=surv[surv.tissue_label_consistent]
assert len(verified)==158 and verified.case_id.is_unique
assert verified.osevent.sum()==56 and verified.rsfevent.sum()==80
score_err=float(np.max(np.abs(surv.protein_OSARS_full-surv.NQO1_self_component-surv.protein_OSARS_without_NQO1)))
assert score_err<1e-10
ps=pd.read_csv(p/'paired_NQO1_statistics.csv').iloc[0]
assert np.isclose(ps.mean_difference,pairs.difference_tumor_minus_adjacent.mean())
assert np.isclose(ps.median_difference,pairs.difference_tumor_minus_adjacent.median())
rho=float(verified.NQO1.rank().corr(verified.protein_OSARS_without_NQO1.rank()))
cs=pd.read_csv(p/'NQO1_surrogate_without_self_correlation.csv').set_index('method')
assert abs(rho-cs.loc['spearman','estimate'])<1e-12
cox=pd.read_csv(p/'NQO1_continuous_cox.csv')
for r in cox.itertuples():
    beta=math.log(r.hr)
    se=(math.log(r.upper)-math.log(r.lower))/(2*1.959963984540054)
    p_normal=math.erfc(abs(beta/se)/math.sqrt(2))
    assert abs(p_normal-r.p)<1e-10
manifest=json.loads((p/'manifest.json').read_text())
for s in manifest['sources']:
    assert hashlib.sha256(Path(s['path']).read_bytes()).hexdigest()==s['sha256']
qa={'status':'PASS','verified_pairs':159,'verified_primary_survival_patients':158,
    'paired_samples_all_originally_measured':True,'tissue_conflicts_preserved_and_excluded':True,
    'surrogate_self_removal_max_abs_error':score_err,'spearman_independently_recalculated':rho,
    'cox_p_values_agree_with_reported_HR_CI':True,'original_source_hashes_unchanged':True,
    'scope':'No imputed pair identities or outcome-selected NQO1 cutoffs'}
(p/'independent_output_qa.json').write_text(json.dumps(qa,indent=2))
print(json.dumps(qa,indent=2))
