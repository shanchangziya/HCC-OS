"""Preserve vector panels at their intended 85 mm width and add letters separately."""
from pathlib import Path
import json
import pymupdf as fitz
from matplotlib import font_manager
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'figures';P=OUT/'panels';MM=72/25.4
font=str(font_manager.findfont(font_manager.FontProperties(family='Arial',weight='bold'),fallback_to_default=True))
d=fitz.open();page=d.new_page(width=170*MM,height=155*MM);page.insert_font(fontname='ArialBold',fontfile=font)
items=[('a',0,0,85,95),('b',85,0,85,95),('c',0,95,85,60),('d',85,95,85,60)]
for letter,x,y,w,h in items:
    src=fitz.open(P/f'sfig_pathology_{letter}.pdf')
    assert abs(src[0].rect.width-w*MM)<.05 and abs(src[0].rect.height-h*MM)<.05
    page.show_pdf_page(fitz.Rect(x*MM,y*MM,(x+w)*MM,(y+h)*MM),src,0)
    page.insert_text(fitz.Point((x+1.2)*MM,(y+4)*MM),letter.upper(),fontname='ArialBold',fontsize=10)
pdf=OUT/'Supplementary_pathology_corrected_internal.pdf';d.save(pdf,garbage=4,deflate=True)
page.get_pixmap(dpi=300,alpha=False).save(pdf.with_suffix('.png'))
spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans']]
outside=[s['text'] for s in spans if not page.rect.contains(fitz.Rect(s['bbox']))]
assert not outside,outside
fonttypes=sorted({f[2] for f in page.get_fonts(full=True)})
assert 'Type3' not in fonttypes and len(page.get_images())==0
qc={'width_mm':170,'height_mm':155,'min_font_pt':min(s['size'] for s in spans),'font_types':fonttypes,
    'raster_images':len(page.get_images()),'text_outside_page':outside,'curve_horizon_years':5,
    'statistics_recomputed_by_plotter':False}
(OUT/'vector_qc.json').write_text(json.dumps(qc,indent=2));print(json.dumps(qc,indent=2))
