"""Recover published-signature genes from unfiltered ICGC FPKM source.

The archived ICGC matrix dropped every gene having any NA across 445 samples.
We retain observed values and explicit NA; no inferred zero or imputation.
"""
from pathlib import Path
import csv,json,os
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
raw=ROOT/'data/raw'
source_value=os.environ.get('HCC_OS_ICGC_FPKM')
if not source_value:
    raise SystemExit('Set HCC_OS_ICGC_FPKM to Merge_RNAseq_FKPM_Symbol.txt.')
source=Path(source_value).expanduser().resolve()
defs=json.loads((ROOT/'references/signature_definitions.json').read_text())
genes=set().union(*(set(d.get('coefficients',{})) for d in defs.values()))
ids=pd.read_csv(raw/'ICGC-LIRI_frozen_risks.csv').ID.tolist()
with source.open() as f:
    reader=csv.reader(f,delimiter='\t');header=next(reader)
    wanted=[header.index(i) for i in ids];rows={}
    for r in reader:
        if r[0] in genes:rows[r[0]]=[float(r[j]) if r[j] not in ('NA','','NaN') else np.nan for j in wanted]
x=pd.DataFrame(rows,index=ids)
assert set(genes)==set(x.columns)
assert not (x<0).any().any()
x=np.log2(x+1)
x.index.name='ID';x.to_csv(raw/'ICGC_unfiltered_published_signature_log_expression.csv')
old=pd.read_csv(raw/'ICGC_full_expression_logscale.csv.gz',index_col=0).T
common=sorted(set(x.columns)&set(old.columns))
error=np.max(np.abs(old.loc[ids,common].to_numpy()-x[common].to_numpy()))
assert error<1e-10
out=pd.DataFrame({'gene':x.columns,'missing_n':x.isna().sum().to_numpy(),'available_n':x.notna().sum().to_numpy()})
out.to_csv(ROOT/'data/processed/ICGC_unfiltered_signature_missingness.csv',index=False)
print('Matched available-gene max absolute difference:',error)
print(out.to_string(index=False))
