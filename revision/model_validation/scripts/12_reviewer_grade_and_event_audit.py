"""Export existing grade-adjusted Cox results and audit competing-event readiness.
No new model fitting; does not alter main analyses or figures.
"""
from pathlib import Path
import hashlib,os
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/processed'
STUDY=Path(os.environ.get('HCC_OS_LEGACY_PROJECT_ROOT',ROOT.parents[1])).expanduser().resolve()

def main():
    cox=pd.read_csv(OUT/'Cox_forest_long.csv')
    ph=pd.read_csv(OUT/'Cox_PH_diagnostics.csv')
    sel=(cox.cohort=='TCGA-LIHC')&(cox.model=='multivariable_plus_grade')
    result=cox.loc[sel].copy()
    diagnostics=ph[(ph.cohort=='TCGA-LIHC')&(ph.model=='multivariable_plus_grade')]
    result=result.merge(diagnostics[['term','p']].rename(columns={'p':'PH_p'}),on='term',validate='one_to_one')
    result['global_PH_p']=float(diagnostics.loc[diagnostics.term=='GLOBAL','p'].iloc[0])
    result['stage_coding']='III/IV versus I/II'
    result['grade_coding']='G3/G4 versus G1/G2'
    result['interpretation']='existing sensitivity, training apparent association; OSARS-specific PH violation; no new model selection'
    result.to_csv(OUT/'reviewer_R2M3_TCGA_grade_sensitivity.csv',index=False)
    d=pd.read_csv(OUT/'TCGA-LIHC_clinical_analysis_rows.csv')
    source=Path(os.environ.get('HCC_OS_TCGA_CLINICAL_GRADE',STUDY/'FIG4_clinical_utility/Fig4_TCGA_OSARS_clinical_all_matched.csv'))
    grade=pd.read_csv(source)[['ID','Grade']].rename(columns={'Grade':'source_Grade'})
    d=d.merge(grade,on='ID',validate='one_to_one')
    assert d.grade.fillna('NA').equals(d.source_Grade.fillna('NA'))
    terms=['time','event','OSARS_SD','age10','male','stage_advanced','grade_high']
    inc=d[['ID','event','grade','stage_num']+terms[2:]].copy()
    inc['included']=d[terms].notna().all(axis=1)
    inc['missing_fields']=d[terms].apply(lambda x:';'.join(x.index[x.isna()]),axis=1)
    assert inc.included.sum()==309 and inc.loc[inc.included,'event'].sum()==104
    inc.to_csv(OUT/'reviewer_R2M3_TCGA_grade_completecase_inclusion.csv',index=False)
    pdc=ROOT/'data/raw/CPTAC_clinical_source.csv'
    x=pd.read_csv(pdc).query("ID!='T724'").copy()
    assert len(x)==158
    x['time_relation']=np.where(np.isclose(x.rsftime,x.ostime),'equal',np.where(x.rsftime<x.ostime,'RFS earlier','RFS later'))
    x.groupby(['osevent','rsfevent','time_relation']).size().reset_index(name='n').to_csv(OUT/'PDC_competing_event_crossclassification.csv',index=False)
    ambiguity=x[(x.osevent.eq(1)&x.rsfevent.eq(0))|(x.rsftime>x.ostime+1e-6)].copy()
    ambiguity['audit_reason']=np.where(ambiguity.rsftime>ambiguity.ostime+1e-6,'RFS time exceeds recorded OS follow-up; endpoint follow-up basis unresolved','death with no recorded RFS event: possible competing event, subject to endpoint-definition confirmation')
    ambiguity[['ID','ostime','osevent','rsftime','rsfevent','time_relation','audit_reason']].to_csv(OUT/'PDC_competing_event_ambiguities.csv',index=False)
    fields=[
      ('OS time/status',True,'ostime;osevent','present; months and binary status as used by archived survival scripts'),
      ('RFS time/status',True,'rsftime;rsfevent','present; exact original event definition not available in the archived table'),
      ('validated three-state recurrence/death/censor code',False,'','absent'),
      ('separate recurrence date and death date',False,'','absent'),
      ('recurrence ascertainment / endpoint data dictionary',False,'','not found in examined local protein sources'),
      ('death without prior recurrence confirmed as competing event',False,'','17 rows have osevent=1 and rsfevent=0; classifying them as competing events requires endpoint-definition confirmation'),
      ('consistent endpoint follow-up basis',False,'','T615: recorded RFS time 38.93 months exceeds recorded OS censoring time 37.37 months; clarify before a joint event process'),
      ('Fine-Gray or cumulative-incidence analysis ready',False,'','not performed; no verified cause-coded event process can be reconstructed without additional assumptions')]
    pd.DataFrame(fields,columns=['required_information','available','fields','finding']).to_csv(OUT/'PDC_competing_event_field_audit.csv',index=False)
    paths=[source,ROOT/'data/raw/TCGA_clinical_source.csv',pdc,OUT/'Cox_forest_long.csv',OUT/'Cox_PH_diagnostics.csv',Path(__file__)]
    rows=[]
    for f in paths:
        try:path_label=str(f.resolve().relative_to(ROOT.parents[1]))
        except ValueError:path_label=f'external:{f.name}'
        rows.append(dict(path=path_label,bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
    pd.DataFrame(rows).to_csv(ROOT/'reviewer_followup_source_manifest.csv',index=False)
    print('Confirmed existing grade+stage sensitivity: n=309, events=104. PDC competing-risk fields audited; no model refitted.')

if __name__=='__main__':main()
