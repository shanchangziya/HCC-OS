"""Reproduce analysis from exported data/raw inputs, using an explicit Rscript.
Run from the manuscript project directory. Source export is documented in README.
"""
from pathlib import Path
import subprocess,sys,os
ROOT=Path(__file__).resolve().parents[1]
R=os.environ.get('HCC_OS_RSCRIPT','Rscript')
for name in ['04_nested_tcga_cv.py','05_prepare_scores_and_metrics.py','06_clinical_survival.R',
             '07_nested_and_incremental_summaries.py','08_sensitivity_statistics.py',
             '09_cox_assumption_sensitivities.R','10_verify_and_manifest.py']:
    exe=R if name.endswith('.R') else sys.executable
    print('Running',name,flush=True)
    subprocess.run([exe,str(ROOT/'scripts'/name)],check=True)
