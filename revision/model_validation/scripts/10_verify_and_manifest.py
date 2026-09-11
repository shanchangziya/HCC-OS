"""Validate statistically material invariants and inventory reproducibility files."""
from pathlib import Path
import json,hashlib,platform,os
import numpy as np
import pandas as pd
from survival_metrics import bootstrap_metrics,cindex

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/processed';RAW=ROOT/'data/raw'

def main():
    checks=[]
    def check(name,value,detail=''):
        assert value,name
        checks.append(dict(check=name,passed=bool(value),detail=detail))
    a=pd.read_csv(OUT/'timeROC_point_estimate_reproduction_check.csv')
    check('all six archived timeROC point estimates reproduced',len(a)==6 and a.absolute_error.max()<1e-12,str(a.absolute_error.max()))
    r=pd.read_csv(OUT/'univ_test_R.csv');p=pd.read_csv(OUT/'univ_test_python.csv')
    check('vectorized Efron Cox equals survival::coxph on five actual genes',np.allclose(r[['beta','p']],p[['beta','p']],atol=1e-12))
    outer=pd.read_csv(OUT/'nested_outer_sample_splits.csv');inner=pd.read_csv(OUT/'nested_inner_sample_splits.csv')
    for fold,d in outer.groupby('outer_fold'):
        tr=set(d.loc[d.role=='train','ID']);te=set(d.loc[d.role=='test','ID'])
        check(f'outer {fold} train and test disjoint',not tr&te and len(tr|te)==343)
        for f,dd in inner[inner.outer_fold==fold].groupby('inner_fold'):
            itr=set(dd.loc[dd.role=='train','ID']);iva=set(dd.loc[dd.role=='validation','ID'])
            check(f'outer {fold} inner {f} no held-out contamination',not itr&iva and itr|iva==tr and not (itr|iva)&te)
    test=outer[outer.role=='test']
    check('each TCGA sample held out exactly once',len(test)==343 and test.ID.is_unique)
    t=pd.read_csv(OUT/'TCGA-LIHC_all_scores.csv')
    s,b,cube=bootstrap_metrics(t.time,t.event,{'risk':t.OSARS_frozen},(),20,20260905)
    rng=np.random.default_rng(20260905);direct=[]
    for _ in range(20):
        ix=rng.integers(len(t),size=len(t));d=t.iloc[ix];direct.append(cindex(d.time,d.event,d.OSARS_frozen))
    check('optimized bootstrap equals explicit row-resampling Harrell C',np.allclose(cube[:,0,0],direct,atol=1e-14))
    c=pd.read_csv(OUT/'ICGC-LIRI_benchmark_completecases_all_scores.csv')
    check('published models compared on same 141 complete samples',len(c)==141 and c[['OSARS_frozen','Hong2022_8gene_log','Ma2024_TR_OSRG_3gene_zscore']].notna().all().all())
    d=pd.read_csv(OUT/'PDC000198_OS_verified158_all_scores.csv')
    check('verified PDC tumor-only primary excludes T724',len(d)==158 and 'T724' not in d.ID.to_list())
    dup=pd.read_csv(OUT/'ICGC_identical_expression_profile_audit.csv')
    cl=pd.read_csv(RAW/'ICGC_clinical_source.csv')
    dup[dup.identical_expression_profile_n>1].merge(cl,on='ID',validate='one_to_one').to_csv(OUT/'ICGC_duplicate_expression_clinical.csv',index=False)
    # ICGC donor IDs are unavailable, hence interval labels consistently use rows.
    for f in OUT.glob('*.csv'):
        text=f.read_text()
        if 'paired patient percentile bootstrap; fixed scores' in text:
            f.write_text(text.replace('paired patient percentile bootstrap; fixed scores','paired sample-row percentile bootstrap; fixed scores'))
    (OUT/'verification_report.json').write_text(json.dumps(checks,indent=2))
    discovery_root=Path(os.environ.get('HCC_OS_DISCOVERY_ROOT',ROOT.parents[1])).expanduser().resolve()
    source_paths=[discovery_root/'QWEN0208_no_GSE14520.Rdata',discovery_root/'res_no_GSE14520.Rdata',
      discovery_root/'StepCox_GBM_model_object.Rdata',
      ROOT.parent/'single_cell/data/processed/discovery_NMF_genes.csv',
      ROOT.parent/'immune_therapy/data/processed/immune_patient_values.csv',
      ROOT.parent/'protein_support/data/processed/pdc_source_tissue_discrepancies.csv']
    paths=source_paths+list(RAW.glob('*'))+list((ROOT/'scripts').glob('*.py'))+list((ROOT/'scripts').glob('*.R'))+list((ROOT/'references').glob('*'))
    manifest=[]
    for path in paths:
        if not path.is_file():continue
        sha=hashlib.sha256()
        with path.open('rb') as handle:
            for block in iter(lambda:handle.read(1024*1024),b''):sha.update(block)
        try:path_label=str(path.resolve().relative_to(ROOT.parents[1]))
        except ValueError:path_label=f'external:{path.name}'
        manifest.append(dict(path=path_label,bytes=path.stat().st_size,sha256=sha.hexdigest()))
    pd.DataFrame(manifest).to_csv(ROOT/'source_and_code_manifest.csv',index=False)
    print(f'{len(checks)} statistical verification checks passed; {len(manifest)} source/code files inventoried.')

if __name__=='__main__':main()
