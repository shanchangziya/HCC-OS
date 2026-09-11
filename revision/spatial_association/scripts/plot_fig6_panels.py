"""Descriptive Fig6 from analyst-verified processed spatial tables.

No correlation, testing, abundance normalization, or image interpolation is
performed here. Map radii are supplied geometric display suggestions.
"""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib import font_manager

ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data/processed';OUT=ROOT/'figures/panels'
OUT.mkdir(parents=True,exist_ok=True);MM=1/25.4
mpl.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':8,
 'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.75,
 'xtick.major.width':.75,'ytick.major.width':.75,'lines.linewidth':1,'pdf.fonttype':42,
 'ps.fonttype':42,'svg.fonttype':'none','figure.dpi':300,'savefig.dpi':300,
 'savefig.bbox':'tight','savefig.pad_inches':.02})
font=str(font_manager.findfont('Arial',fallback_to_default=True));sources={};records=[]
def read(name):
 p=DATA/name;sources[name]={'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()};return pd.read_csv(p)
def export(fig,name,w,h,inputs):
 with mpl.rc_context({'savefig.bbox':None}):
  for ext in ['pdf','png']:fig.savefig(OUT/(name+'.'+ext),dpi=300,bbox_inches=None)
 records.append({'name':name,'width_mm':w,'height_mm':h,'inputs':inputs,'font':font,'min_font_pt':7});plt.close(fig)
data=read('spatial_spot_coordinates_and_abundance.csv');meta=read('section_metadata_and_source.csv');rho=read('section_TAM_descriptive_spearman.csv');display=read('pooled_display_quantiles.csv')
sections=['HCC-1T','HCC-2T','HCC-3T','HCC-4T'];tams=['TREM2 TAM','SLC40A1 TAM','HSP TAM']

# A: all twelve supplied within-section correlations, without inferential marks.
matrix=rho.pivot(index='tam_label',columns='section',values='spearman_rho').loc[tams,sections]
fig=plt.figure(figsize=(170*MM,50*MM));ax=fig.add_axes([.21,.29,.65,.58]);cax=fig.add_axes([.90,.29,.018,.58])
fig.text(.535,.948,'Historical OS-high–TAM abundance',ha='center',fontsize=8)
mesh=ax.pcolormesh(np.arange(5)-.5,np.arange(4)-.5,matrix.to_numpy(),cmap='RdBu_r',vmin=-1,vmax=1,edgecolors='white',linewidth=.8,rasterized=False)
for i in range(3):
 for j in range(4):
  v=matrix.iloc[i,j];ax.text(j,i,f'{v:.3f}',ha='center',va='center',color='white' if abs(v)>.6 else '#222222',fontsize=7)
ax.set_yticks(range(3),tams);ax.set_xticks(range(4),[f'{s}\n{int(meta[meta.section==s].iloc[0].n_spots):,} spots' for s in sections]);ax.invert_yaxis();ax.tick_params(length=0,pad=5)
for sp in ax.spines.values():sp.set_visible(False)
cb=fig.colorbar(mesh,cax=cax,ticks=[-1,0,1]);cb.set_label('Spearman ρ',fontsize=8,labelpad=3);cb.outline.set_visible(False);cb.ax.tick_params(length=2,pad=2)
cb.solids.set_rasterized(False)
export(fig,'fig6a',170,50,['section_TAM_descriptive_spearman.csv','section_metadata_and_source.csv'])

# B-E: one identical abundance color scale for every displayed state and section.
# Fixed limits 0..0.5 encompass all exported values without clipping.
assert (display.minimum>=0).all() and (display.maximum<=.5).all()
norm=mpl.colors.Normalize(vmin=0,vmax=.5)
for letter,section in zip('bcde',sections):
 d=data[data.section==section];m=meta[meta.section==section].iloc[0]
 assert len(d)==int(m.n_spots)
 fig=plt.figure(figsize=(85*MM,75*MM));fig.text(.5,.952,section,ha='center',fontsize=8)
 for j,(state,label) in enumerate([('Tumor_OS-high','Historical OS-high'),('Mac_TREM2_TAM','TREM2 TAM')]):
  ax=fig.add_axes([.04+j*.48,.23,.44,.60]);r=m.map_point_radius_suggestion
  circles=[Circle((x,y),radius=r) for x,y in zip(d.imagecol,d.imagerow)]
  pc=PatchCollection(circles,cmap='viridis',norm=norm,edgecolor='none',linewidth=0,rasterized=False)
  pc.set_array(d[state].to_numpy());ax.add_collection(pc)
  ax.set_xlim(m.x_min-1.5*r,m.x_max+1.5*r);ax.set_ylim(m.y_max+1.5*r,m.y_min-1.5*r)
  ax.set_aspect('equal',adjustable='box');ax.axis('off');fig.text(.26+j*.48,.862,label,ha='center',fontsize=7)
 cax=fig.add_axes([.22,.14,.56,.035]);cb=fig.colorbar(pc,cax=cax,orientation='horizontal',ticks=[0,.25,.5])
 cb.set_label('Posterior mean abundance',fontsize=8,labelpad=3);cb.ax.tick_params(length=2,pad=2);cb.outline.set_visible(False)
 cb.solids.set_rasterized(False)
 export(fig,'fig6'+letter,85,75,['spatial_spot_coordinates_and_abundance.csv','section_metadata_and_source.csv','pooled_display_quantiles.csv'])

(ROOT/'figures/panel_manifest.json').write_text(json.dumps({'panels':records,'sources':sources,'color_limits':[0,.5],'abundance_clipped':False,'display_radius_source':'section_metadata_and_source.csv: map_point_radius_suggestion','spatial_interpolation':False},indent=2))
print('Exported five Fig6 panels from all four sections.')
