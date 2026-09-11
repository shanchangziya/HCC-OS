"""Compose new Fig6 vector panels at their intended physical sizes."""
from pathlib import Path
import json
import pymupdf as fitz
from matplotlib import font_manager
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'figures';MM=72/25.4
font=str(font_manager.findfont(font_manager.FontProperties(family='Arial',weight='bold'),fallback_to_default=True))
doc=fitz.open();p=doc.new_page(width=170*MM,height=200*MM);p.insert_font(fontname='ArialBold',fontfile=font)
for letter,x,y,w,h in [('a',0,0,170,50),('b',0,50,85,75),('c',85,50,85,75),('d',0,125,85,75),('e',85,125,85,75)]:
 src=fitz.open(OUT/'panels'/('fig6'+letter+'.pdf'))
 assert abs(src[0].rect.width-w*MM)<.05 and abs(src[0].rect.height-h*MM)<.05
 p.show_pdf_page(fitz.Rect(x*MM,y*MM,(x+w)*MM,(y+h)*MM),src,0)
 p.insert_text(fitz.Point((x+1.2)*MM,(y+4)*MM),letter.upper(),fontname='ArialBold',fontsize=10)
f=OUT/'Fig6_descriptive_spatial_associations.pdf';doc.save(f,garbage=4,deflate=True);p.get_pixmap(dpi=300,alpha=False).save(f.with_suffix('.png'))
spans=[s for b in p.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans']]
outside=[s['text'] for s in spans if not p.rect.contains(fitz.Rect(s['bbox']))]
qc={'width_mm':170,'height_mm':200,'minimum_font_pt':min(s['size'] for s in spans),'font_types':sorted({f[2] for f in p.get_fonts(full=True)}),'raster_images':len(p.get_images()),'outside_page':outside}
(OUT/'vector_qc.json').write_text(json.dumps(qc,indent=2));print(json.dumps(qc,indent=2))
assert not outside and qc['minimum_font_pt']>=7 and qc['raster_images']==0
