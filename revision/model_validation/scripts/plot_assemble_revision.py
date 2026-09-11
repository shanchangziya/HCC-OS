"""Compose newly rendered, editable vector panels without rescaling them."""
from pathlib import Path
import json
import pymupdf as fitz
from matplotlib import font_manager
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'figures';MM=72/25.4
font=str(font_manager.findfont(font_manager.FontProperties(family='Arial',weight='bold'),fallback_to_default=True))
specs={
 'Fig3_program_specificity_and_model_validation':(170,229,[
 ('fig3a','A',0,0,85,62),('fig3b','B',85,0,85,62),
 ('fig3c','C',0,62,85,80),('fig3d','D',85,62,85,80),
 ('fig3e','E',0,142,85,87),('fig3f','F',85,142,85,87)]),
 'Fig4_survival_and_clinical_association':(170,163,[
 ('fig4a','A',0,0,170/3,95),('fig4b','B',170/3,0,170/3,95),('fig4c','C',340/3,0,170/3,95),
 ('fig4d','D',0,95,85,68),('fig4e','E',85,95,85,68)]),
 'Supplementary_published_signature_AUC':(170,70,[('sfig_auc_a','A',0,0,85,70),('sfig_auc_b','B',85,0,85,70)]),
 'Supplementary_clinical_sensitivity':(170,160,[
 ('sfig_clinical_a','A',0,0,85,95),('sfig_clinical_b','B',85,0,85,95),
 ('sfig_clinical_c','C',0,95,85,65),('sfig_clinical_d','D',85,95,85,65)]),
 'Supplementary_Cox_nonlinearity':(170,80,[
 ('sfig_spline_a','A',0,0,170/3,80),('sfig_spline_b','B',170/3,0,170/3,80),('sfig_spline_c','C',340/3,0,170/3,80)])
}
report={}
for name,(w,h,items) in specs.items():
    doc=fitz.open();p=doc.new_page(width=w*MM,height=h*MM);p.insert_font(fontname='ArialBold',fontfile=font)
    for srcname,letter,x,y,pw,ph in items:
        src=fitz.open(OUT/'panels'/(srcname+'.pdf'))
        assert abs(src[0].rect.width-pw*MM)<.05 and abs(src[0].rect.height-ph*MM)<.05
        p.show_pdf_page(fitz.Rect(x*MM,y*MM,(x+pw)*MM,(y+ph)*MM),src,0)
        p.insert_text(fitz.Point((x+1.2)*MM,(y+4)*MM),letter,fontname='ArialBold',fontsize=10)
    path=OUT/(name+'.pdf');doc.save(path,garbage=4,deflate=True);p.get_pixmap(dpi=300,alpha=False).save(path.with_suffix('.png'))
    spans=[s for b in p.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans']]
    outside=[s['text'] for s in spans if not p.rect.contains(fitz.Rect(s['bbox']))]
    qc={'width_mm':w,'height_mm':h,'minimum_font_pt':min(s['size'] for s in spans),'font_types':sorted({f[2] for f in p.get_fonts(full=True)}),'raster_images':len(p.get_images()),'outside_page':outside}
    report[name]=qc
(OUT/'vector_qc.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
assert all(not q['outside_page'] and q['raster_images']==0 and q['minimum_font_pt']>=7 for q in report.values())
