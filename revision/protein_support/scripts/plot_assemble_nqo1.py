"""Compose fresh NQO1 vector panels at 1:1 scale, with lettering added here."""
from pathlib import Path
import json
import pymupdf as fitz
from matplotlib import font_manager
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'figures';MM=72/25.4
font=str(font_manager.findfont(font_manager.FontProperties(family='Arial',weight='bold'),fallback_to_default=True))
doc=fitz.open();p=doc.new_page(width=170*MM,height=85*MM);p.insert_font(fontname='ArialBold',fontfile=font)
for i,letter in enumerate('abc'):
    src=fitz.open(OUT/'panels'/f'sfig_nqo1_{letter}.pdf');w=170/3
    assert abs(src[0].rect.width-w*MM)<.05
    p.show_pdf_page(fitz.Rect(i*w*MM,0,(i+1)*w*MM,85*MM),src,0)
    p.insert_text(fitz.Point((i*w+1.2)*MM,4*MM),letter.upper(),fontname='ArialBold',fontsize=10)
f=OUT/'Supplementary_NQO1_orthogonal_protein.pdf';doc.save(f,garbage=4,deflate=True);p.get_pixmap(dpi=300,alpha=False).save(f.with_suffix('.png'))
spans=[s for b in p.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans']]
out=[s['text'] for s in spans if not p.rect.contains(fitz.Rect(s['bbox']))]
assert not out,out
qc={'width_mm':170,'height_mm':85,'minimum_font_pt':min(s['size'] for s in spans),'font_types':sorted({f[2] for f in p.get_fonts(full=True)}),'raster_images':len(p.get_images()),'outside_page':out}
assert qc['minimum_font_pt']>=7 and qc['raster_images']==0
(OUT/'vector_qc.json').write_text(json.dumps(qc,indent=2));print(json.dumps(qc,indent=2))
