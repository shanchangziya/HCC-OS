"""Verify 67 eligible BIOSTORM labels against downloaded original GEO annotations."""
from pathlib import Path
import hashlib,json,re
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
text=(ROOT/'data/source/GSE109211_sample_annotations.txt').read_text()
rows=[]
for block in text.split('^SAMPLE = ')[1:]:
    sid=block.splitlines()[0].strip()
    treatment=re.search(r'!Sample_characteristics_ch1 = treatment: ([^\n\r]+)',block)
    outcome=re.search(r'!Sample_characteristics_ch1 = outcome: ([^\n\r]+)',block)
    if treatment and outcome:rows.append(dict(ID=sid,GEO_treatment=treatment[1],GEO_outcome=outcome[1]))
geo=pd.DataFrame(rows)
geo.to_csv(ROOT/'data/processed/GSE109211_original_GEO_labels.csv',index=False)
clinical=pd.read_csv(ROOT/'data/processed/therapy_patients.csv')
merged=clinical[clinical.Cohort.str.contains('GSE109211')].merge(geo,on='ID',validate='one_to_one')
assert len(merged)==67 and (merged.GEO_treatment=='Sor').all()
assert (merged.Response==merged.GEO_outcome.map({'responder':'R','non-responder':'NR'})).all()
summary='All 67 eligible sorafenib cases match original GEO sample ID, treatment=Sor, and outcome labels exactly. Original BIOSTORM publication identifies these as recurrence-free-survival-benefit signature classes, not radiographic tumor response.\n'
(ROOT/'qc/GEO_label_validation.txt').write_text(summary)
manifest=[]
for f in (ROOT/'data/source').iterdir():
    if f.suffix in ['.txt','.xml']:manifest.append({'file':f.name,'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'retrieved':'2026-09-05'})
(ROOT/'qc/public_source_manifest.json').write_text(json.dumps(manifest,indent=2))
print(summary)
