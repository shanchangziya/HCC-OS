"""Separate typesetting step: compose fresh vector PDFs at 1:1 physical scale."""
from pathlib import Path
import json
import pymupdf as fitz
from matplotlib import font_manager
ROOT=Path(__file__).resolve().parents[1];PANELS=ROOT/'figures/panels';OUT=ROOT/'figures/composites'
MM=72/25.4
fontfile=str(font_manager.findfont(font_manager.FontProperties(family='Arial',weight='bold'),fallback_to_default=True))
qc=[]

def compose(name,width,height,items):
    doc=fitz.open();page=doc.new_page(width=width*MM,height=height*MM)
    page.insert_font(fontname='ArialBold',fontfile=fontfile)
    for stem,x,y,w,h,letter in items:
        src=fitz.open(PANELS/(stem+'.pdf'))
        assert abs(src[0].rect.width-w*MM)<.05 and abs(src[0].rect.height-h*MM)<.05
        page.show_pdf_page(fitz.Rect(x*MM,y*MM,(x+w)*MM,(y+h)*MM),src,0)
        if letter:page.insert_text(fitz.Point((x+1.2)*MM,(y+4)*MM),letter,fontname='ArialBold',fontsize=10)
    doc.save(OUT/(name+'.pdf'),garbage=4,deflate=True)
    page.get_pixmap(dpi=300,alpha=False).save(OUT/(name+'.png'))
    fonts=page.get_fonts(full=True);types=sorted({f[2] for f in fonts})
    images=page.get_images(full=True)
    assert 'Type3' not in types and len(images)==0
    text=page.get_text('dict');spans=[s for b in text['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans']]
    outside=[s['text'] for s in spans if not page.rect.contains(fitz.Rect(s['bbox']))]
    assert not outside, outside
    qc.append(dict(file=name,width_mm=width,height_mm=height,font_types=types,embedded_fonts=len(fonts),raster_images=len(images),min_font_pt=min(s['size'] for s in spans),text_outside_page=outside))
    doc.close()

OUT.mkdir(exist_ok=True,parents=True)
compose('Fig5_immune_TIDE_revised',170,125,[('fig5a',0,0,85,77,'A'),('fig5b',85,0,85,125,'B'),('fig5c',0,77,85,48,'C')])
compose('Supplementary_therapy_contexts',170,68,[(f'sfig_therapy_{l}',i*170/3,0,170/3,68,l.upper()) for i,l in enumerate('abc')])
compose('Supplementary_DrugReflector',85,110,[('sfig_drug_a',0,0,85,110,'A')])
compose('Supplementary_checkpoints',170,95,[('sfig_checkpoints_a',0,0,170,95,'A')])
compose('Supplementary_CIBERSORT_fit_sensitivity',170,63,[('sfig_cibersort_a',0,0,85,63,'A'),('sfig_cibersort_b',85,0,85,63,'B')])
(ROOT/'qc/vector_export_qc.json').write_text(json.dumps(qc,indent=2))
print(json.dumps(qc,indent=2))
