#!/usr/bin/env python3
"""Prepare verified PDC paired NQO1 data and numerical survival analyses; no plots."""
import argparse, csv, hashlib, json, os, subprocess, urllib.request
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def digest(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--proteomics-root',default=os.environ.get('HCC_OS_PROTEOMICS_ROOT'))
    ap.add_argument('--surrogate-weights',default=os.environ.get('HCC_OS_PDC_SURROGATE_WEIGHTS'))
    ap.add_argument('--surrogate-scores',default=os.environ.get('HCC_OS_PDC_SURROGATE_SCORES'))
    ap.add_argument('--rscript',default=os.environ.get('HCC_OS_RSCRIPT','Rscript'))
    a=ap.parse_args()
    if not a.proteomics_root:
        ap.error('provide --proteomics-root or HCC_OS_PROTEOMICS_ROOT')
    base=Path(a.proteomics_root).expanduser().resolve()
    weights=Path(a.surrogate_weights).expanduser().resolve() if a.surrogate_weights else ROOT/'data/raw/CPTAC_OSARS_protein_surrogate_weights.csv'
    scores=Path(a.surrogate_scores).expanduser().resolve() if a.surrogate_scores else ROOT/'data/raw/CPTAC_OSARS_protein_surrogate_scores_OS_RFS.csv'
    raw=ROOT/'data/raw';processed=ROOT/'data/processed';logs=ROOT/'logs'
    for p in (raw,processed,logs):p.mkdir(parents=True,exist_ok=True)
    sources=[base/'肿瘤和临床数据/CPTACLIHCtumor.Rdata', base/'肝癌和正常组织/CPTACliver.Rdata',
             base/'肝癌和正常组织/CPTACLIVER.txt',base/'肝癌和正常组织/CPTAC04.normalize.R',
             weights,scores]
    missing=[str(p) for p in sources if not p.is_file()]
    if missing:ap.error('missing required source file(s): '+', '.join(missing))
    manifest={'date_utc':datetime.now(timezone.utc).isoformat(),'sources':[{'source_id':p.name,'sha256':digest(p),'bytes':p.stat().st_size} for p in sources]}
    (processed/'manifest.json').write_text(json.dumps(manifest,indent=2))
    query='{biospecimenPerStudy(pdc_study_id:"PDC000198"){case_id case_submitter_id sample_id sample_submitter_id sample_type aliquot_id aliquot_submitter_id} study(pdc_study_id:"PDC000198"){study_name study_description pdc_study_id study_id}}'
    (raw/'pdc_pairing_query.txt').write_text(query)
    cache=raw/'pdc_000198_biospecimen.json'
    if not cache.exists():
        req=urllib.request.Request('https://pdc.cancer.gov/graphql',data=json.dumps({'query':query}).encode(),headers={'Content-Type':'application/json'})
        cache.write_text(json.dumps(json.load(urllib.request.urlopen(req,timeout=30)),indent=2))
    biospecimen=pd.DataFrame(json.loads(cache.read_text())['data']['biospecimenPerStudy'])
    assert biospecimen.aliquot_submitter_id.is_unique
    biospecimen.to_csv(processed/'pdc_biospecimen_mapping.csv',index=False)
    # Export exact source values with R, including the original protein surrogate.
    cmd=[a.rscript,str(ROOT/'scripts/export_protein.R'),str(sources[0]),str(sources[1]),str(sources[4]),str(sources[5]),str(raw)]
    with (logs/'export.log').open('w') as f:subprocess.run(cmd,check=True,stdout=f,stderr=subprocess.STDOUT)
    expression=pd.read_csv(raw/'normalized_nqo1_all_samples.csv',dtype={'sample_id':str})
    expression['aliquot_submitter_id']=expression.sample_id.str.replace(r'^[TP]','',regex=True)
    expression=expression.merge(biospecimen,on='aliquot_submitter_id',how='left',validate='one_to_one',suffixes=('','_pdc'))
    assert expression.case_id.notna().all(), 'Some expression samples cannot be matched to PDC aliquots'
    expression['tissue_label_consistent']=(expression.sample_id.str.startswith('T'))==(expression.sample_type=='Primary Tumor')
    expression.loc[~expression.tissue_label_consistent].to_csv(processed/'pdc_source_tissue_discrepancies.csv',index=False)
    # Recover which normalized values had originally been missing, before impute.knn.
    with sources[2].open() as f:
        reader=csv.reader(f,delimiter='\t');header=next(reader)
        nqo1=next(row for row in reader if row[0]=='NQO1')
    original={h.split(' ')[0]:v for h,v in zip(header,nqo1) if h.endswith('Unshared Log Ratio')}
    expression['raw_unshared_log_ratio']=pd.to_numeric(expression.sample_id.map(original),errors='coerce')
    expression['originally_measured']=expression.raw_unshared_log_ratio.notna()
    assert expression.normalized_NQO1.notna().all()
    expression['tissue']=np.where(expression.sample_id.str.startswith('T'),'Tumor','Adjacent liver')
    expression.to_csv(processed/'nqo1_all_samples_with_verified_mapping.csv',index=False)
    eligible=expression[expression.tissue_label_consistent]
    pairing_counts=eligible.groupby('case_id').agg(n_samples=('sample_id','size'),n_tissues=('tissue','nunique'))
    true_pair_ids=pairing_counts.index[(pairing_counts.n_samples==2)&(pairing_counts.n_tissues==2)]
    paired=expression[expression.case_id.isin(true_pair_ids)]
    wide=paired.pivot(index='case_id',columns='tissue',values='normalized_NQO1').reset_index()
    wide=wide.rename(columns={'Tumor':'tumor','Adjacent liver':'adjacent'})
    wide['difference_tumor_minus_adjacent']=wide.tumor-wide.adjacent
    measurement=paired.groupby('case_id').originally_measured.all()
    wide['both_originally_measured']=wide.case_id.map(measurement)
    wide.to_csv(processed/'nqo1_verified_pairs.csv',index=False)
    expression[~expression.case_id.isin(true_pair_ids)].to_csv(processed/'nqo1_unpaired_excluded_from_paired_test.csv',index=False)
    survival=pd.read_csv(raw/'nqo1_survival_surrogate.csv')
    survival=survival.merge(expression[['sample_id','case_id','originally_measured','tissue_label_consistent']],on='sample_id',validate='one_to_one')
    assert survival.case_id.is_unique
    survival.to_csv(processed/'nqo1_survival_input.csv',index=False)
    qa={'n_normalized_samples':len(expression),'n_tumor':int((expression.tissue=='Tumor').sum()),
        'n_adjacent':int((expression.tissue=='Adjacent liver').sum()),'n_verified_patient_pairs':len(wide),
        'n_pair_excluded_samples':int((~expression.case_id.isin(true_pair_ids)).sum()),
        'n_pairs_both_originally_measured':int(wide.both_originally_measured.sum()),
        'n_tumor_survival':len(survival),'n_os_events':int(survival.osevent.sum()),'n_rfs_events':int(survival.rsfevent.sum()),
        'n_survival_tissue_label_verified':int(survival.tissue_label_consistent.sum()),
        'tissue_label_conflict_samples':expression.loc[~expression.tissue_label_consistent,'sample_id'].tolist(),
        'pairing_rule':'Same PDC case_id with exactly one Primary Tumor and one Solid Tissue Normal aliquot',
        'survival_time_unit':'Original ostime/rsftime scale; unit not explicitly encoded in RData; no conversion; Cox HR invariant to common time scaling',
        'expression_unit':'Unshared log2 ratio after original KNN imputation and quantile normalization',
        'local_cliincal_txt_rejected':'ICGC SA sample IDs do not match this proteomics cohort',
        'surrogate_self_contribution':'NQO1 z-score times its weight removed exactly before correlation',
        'pdc_mapping_sha256':digest(cache)}
    (processed/'input_qa.json').write_text(json.dumps(qa,indent=2))
    with (logs/'statistics.log').open('w') as f:subprocess.run([a.rscript,str(ROOT/'scripts/protein_statistics.R'),str(processed)],check=True,stdout=f,stderr=subprocess.STDOUT)
    manifest['sources_unchanged']=all(digest(p)==item['sha256'] for p,item in zip(sources,manifest['sources']))
    manifest['completed_utc']=datetime.now(timezone.utc).isoformat()
    (processed/'manifest.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps(qa,indent=2))


if __name__=='__main__':main()
