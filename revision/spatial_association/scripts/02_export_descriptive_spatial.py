#!/usr/bin/env python3
"""Export frozen cell2location results; no model refit or inferential spot tests.

Input: ../data/raw copied from the original server cell2location_misty tree.
Output: ../data/processed, with read-only input integrity and numerical audit.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

BASE = Path(__file__).resolve().parents[1]
RAW, OUT = BASE / 'data/raw', BASE / 'data/processed'
OUT.mkdir(parents=True, exist_ok=True)
SECTIONS = ['HCC-1T', 'HCC-2T', 'HCC-3T', 'HCC-4T']
TAM = ['Mac_TREM2_TAM', 'Mac_SLC40A1_TAM', 'Mac_HSP_TAM']
FEATURES = ['Tumor_OS-high', *TAM]
UNIT = 'cell2location posterior mean cell abundance (means_cell_abundance_w_sf)'
LABEL = 'Historical OS_group_valley-derived reference; not revision top-20% labels'

manifest = json.loads((BASE/'audit/source_manifest.json').read_text())
for item in manifest:
    p = RAW / item['relative_path']
    assert p.exists() and hashlib.sha256(p.read_bytes()).hexdigest() == item['sha256'], p

original = pd.read_csv(RAW/'results/colocalization/cell2location_spot_correlations.csv')
run = pd.read_csv(RAW/'results/cell2location_spatial/spatial_run_summary.csv').set_index('sample')
spots, section_rows, correlation_rows, control_rows, qa_rows = [], [], [], [], []
for section in SECTIONS:
    abundance_rel = f'results/cell2location_spatial/{section}/cell_abundance_mean.csv'
    coords_rel = f'inputs/spatial/{section}/spatial_coords.tsv.gz'
    a = pd.read_csv(RAW/abundance_rel, index_col=0)
    c = pd.read_csv(RAW/coords_rel, sep='\t').set_index('barcode')
    assert a.index.is_unique and c.index.is_unique
    assert set(a.index) == set(c.index), section
    assert len(a) == int(run.loc[section, 'n_spots'])
    assert np.isfinite(a.to_numpy()).all() and (a.to_numpy() >= 0).all()
    c = c.loc[a.index]
    assert np.isfinite(c[['imagecol', 'imagerow']].to_numpy()).all()
    assert not c[['imagecol', 'imagerow']].duplicated().any()
    z = c[['imagecol', 'imagerow', 'image']].copy()
    z.insert(0, 'section', section)
    z.insert(1, 'barcode', a.index)
    for f in FEATURES + ['Tumor_OS-low']:
        z[f] = a[f]
    z['total_abundance_all_28_states'] = a.sum(axis=1)
    z.reset_index(drop=True, inplace=True)
    spots.append(z)
    xy = c[['imagecol', 'imagerow']].to_numpy()
    spacing = float(np.median(cKDTree(xy).query(xy, k=2)[0][:, 1]))
    section_rows.append(dict(section=section, n_spots=len(a), n_reference_states=a.shape[1],
        n_model_genes=int(run.loc[section,'n_genes']), max_epochs=int(run.loc[section,'max_epochs']),
        abundance_source=abundance_rel, coordinate_source=coords_rel,
        abundance_unit=UNIT, reference_label_definition=LABEL,
        coordinate_x='imagecol',coordinate_y='imagerow', y_increases='downward',
        coordinate_unit='Seurat exported image coordinates; physical micrometre scale unverified',
        median_nearest_spot_spacing=spacing, map_point_radius_suggestion=0.42*spacing,
        tissue_mask='Observed exported tissue spots only; no interpolated mask',
        x_min=xy[:,0].min(), x_max=xy[:,0].max(), y_min=xy[:,1].min(), y_max=xy[:,1].max(),
        independent_patient_identity='Not established by these exported files'))
    for tumor, mac in [('Tumor_OS-high', m) for m in TAM] + [('Tumor_OS-low','Mac_TREM2_TAM')]:
        rho = float(a[tumor].rank(method='average').corr(a[mac].rank(method='average')))
        saved = original.loc[(original['sample']==section)&(original.tumor==tumor)&(original.mac==mac)]
        assert len(saved)==1 and int(saved.iloc[0].n_spots)==len(a)
        old = float(saved.iloc[0].spearman_rho)
        assert abs(rho-old)<1e-12, (section,tumor,mac,rho,old)
        row = dict(section=section,tumor_state=tumor,macrophage_state=mac,
            tam_label=mac.replace('Mac_','').replace('_TAM',' TAM'), n_spots=len(a),
            spearman_rho=old, statistic='Within-section descriptive Spearman rank correlation',
            inference='No spot-level P value, CI or across-section hypothesis test',
            abundance_unit=UNIT,reference_label_definition=LABEL)
        (correlation_rows if tumor=='Tumor_OS-high' else control_rows).append(row)
        qa_rows.append(dict(section=section,tumor_state=tumor,macrophage_state=mac,
            absolute_rho_difference=abs(rho-old)))

wide = pd.concat(spots,ignore_index=True)
assert len(wide)==16535 and not wide.barcode.duplicated().any()
wide.to_csv(OUT/'spatial_spot_coordinates_and_abundance.csv',index=False)
long = wide.melt(id_vars=['section','barcode','imagecol','imagerow'],value_vars=FEATURES,
                 var_name='cell_state',value_name='posterior_mean_abundance')
long.to_csv(OUT/'spatial_abundance_long.csv',index=False)
sections=pd.DataFrame(section_rows)
sections.to_csv(OUT/'section_metadata_and_source.csv',index=False)
cor=pd.DataFrame(correlation_rows)
cor.to_csv(OUT/'section_TAM_descriptive_spearman.csv',index=False)
cor.pivot(index='section',columns='tam_label',values='spearman_rho').reindex(SECTIONS).to_csv(
    OUT/'section_TAM_descriptive_spearman_matrix.csv')
pd.DataFrame(control_rows).to_csv(OUT/'section_OSlow_TREM2_descriptive_control.csv',index=False)
pd.DataFrame(qa_rows).to_csv(BASE/'audit/correlation_reproduction_check.csv',index=False)
limits=[]
for f in FEATURES:
    values=wide[f]
    limits.append(dict(cell_state=f,minimum=values.min(),maximum=values.max(),
        p01=values.quantile(.01),median=values.median(),p99=values.quantile(.99),
        unit=UNIT, scope='All four sections; suggested consistent per-state display limits only'))
pd.DataFrame(limits).to_csv(OUT/'pooled_display_quantiles.csv',index=False)
checks=dict(status='PASS',n_sections=4,n_spots=len(wide),n_primary_descriptive_correlations=12,
    max_rho_reproduction_error=max(x['absolute_rho_difference'] for x in qa_rows),
    original_abundance_values_unchanged=True,unique_coordinates_and_exact_barcode_alignment=True,
    no_inferential_spot_P_values_or_CIs_exported=True,no_refit_or_label_redefinition=True,
    original_source_hashes_unchanged=True)
(BASE/'audit/export_qa.json').write_text(json.dumps(checks,indent=2))
print(json.dumps(checks,indent=2))
