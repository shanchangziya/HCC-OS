"""Discovery figure reconstructed from saved cell scores and audited summaries."""
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[1]; P=ROOT/'data/processed'; OUT=ROOT/'figures'
mpl.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':7,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.65,'pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':300})
LOW='#0072B2';HIGH='#D55E00';GRAY='#777777'
D=pd.read_csv(P/'discovery_validation_ready.csv')
S=pd.read_csv(P/'discovery_historical_and_top20_counts.csv')
H=pd.read_csv(P/'discovery_score_histogram.csv')
A=pd.read_csv(P/'discovery_definition_audit.csv').set_index('metric').value
OVERVIEW=[(.92,'GSE202642: 10,413 epithelial discovery cells from 11 tissue samples',8),(.65,'7 tumor and 4 adjacent tissue samples; saved labels retained for original model provenance',7),(.36,'OS composite = [z(AddModuleScore) + z(UCell) + z(AUCell)] / 3',8),(.13,'Historical OS-high: 2,437 cells     |     Revised cohort-wide top 20%: 2,083 cells',7)]

def clean(ax):
    ax.spines[['top','right']].set_visible(False);ax.tick_params(length=2.5,pad=2)
def mapbase(ax):
    ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([])
    for sp in ax.spines.values():sp.set_visible(False)
    ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2')
def score(ax):
    x=D.sort_values('OS_composite')
    q=ax.scatter(x.UMAP1,x.UMAP2,c=x.OS_composite,cmap='viridis',s=.6,linewidths=0,rasterized=True)
    cb=ax.figure.colorbar(q,ax=ax,fraction=.045,pad=.025,shrink=.7);cb.set_label('OS composite');cb.ax.tick_params(length=2)
    mapbase(ax)
def groups(ax):
    for label,c in [('OS-low',LOW),('OS-high',HIGH)]:
        x=D[D.OS_group_valley==label];ax.scatter(x.UMAP1,x.UMAP2,s=.6,linewidths=0,color=c,rasterized=True)
    mapbase(ax)
    hs=[Line2D([],[],ls='',marker='o',ms=3,color=col,label=f'{label} (n = {int((D.OS_group_valley==label).sum()):,})') for label,col in [('OS-low',LOW),('OS-high',HIGH)]]
    ax.legend(handles=hs,loc='upper center',bbox_to_anchor=(.5,-.10),frameon=False,ncol=1)
def histogram(ax):
    for label,c in [('OS-low',LOW),('OS-high',HIGH)]:
        h=H[H.group==label]
        ax.bar(h.bin_left,h.n_cells,width=h.bin_right-h.bin_left,align='edge',color=c,alpha=.85,lw=0)
    ax.axvline(A['top20_cutoff'],color='#222222',lw=.9,ls='--')
    ax.text(.97,.96,'Revised 80th percentile\n1.0173; n high = 2,083',ha='right',va='top',transform=ax.transAxes,fontsize=7)
    ax.set_xlabel('OS composite');ax.set_ylabel('Number of cells');clean(ax)
def fractions(ax):
    x=np.arange(len(S));w=.36
    ax.bar(x-w/2,S.historical_fraction,w,color=HIGH,label='Historical')
    ax.bar(x+w/2,S.top20_fraction,w,color=LOW,label='Revised top 20%')
    ax.set_xticks(x,S['sample'],rotation=45);ax.set_ylim(0,.58);ax.set_ylabel('OS-high fraction');ax.set_xlabel('Tissue sample')
    ax.legend(loc='upper right',ncol=2,frameon=False,handlelength=1,columnspacing=.8);clean(ax)
def agreement(ax):
    x=pd.read_csv(P/'discovery_group_agreement.csv',index_col=0)
    ax.imshow(x,cmap='Blues',aspect='equal')
    for i in range(2):
        for j in range(2):ax.text(j,i,f'{x.iloc[i,j]:,}',ha='center',va='center',color='white' if x.iloc[i,j]>4000 else '#111111')
    ax.set_xticks([0,1],['OS-low','OS-high']);ax.set_yticks([0,1],['OS-low','OS-high'])
    ax.set_xlabel('Revised top-20% definition');ax.set_ylabel('Historical saved definition')
    for s in ax.spines.values():s.set_visible(False)
def save(fig,name,folder):
    fig.savefig(OUT/folder/(name+'.pdf'),facecolor='white');fig.savefig(OUT/folder/(name+'.png'),facecolor='white',dpi=300);plt.close(fig)
funcs=[score,groups,histogram,fractions,agreement]
for i,fn in enumerate(funcs):
    fig,ax=plt.subplots(figsize=(3.35,2.75));fig.subplots_adjust(left=.21,right=.86,bottom=.25,top=.85);fn(ax);save(fig,'Fig1'+chr(66+i),'panels')
fig=plt.figure(figsize=(170/25.4,30/25.4));ax=fig.add_axes([.10,.08,.83,.85]);ax.axis('off')
for y,label,size in OVERVIEW:ax.text(0,y,label,fontsize=size,va='top')
save(fig,'Fig1A','panels')
fig=plt.figure(figsize=(170/25.4,190/25.4))
# Panel A makes score definition and analytical provenance explicit.
ax=fig.add_axes([.10,.865,.83,.105]);ax.axis('off')
ax.text(0,1.06,'A',fontsize=10,fontweight='bold',va='bottom')
for y,label,size in OVERVIEW:ax.text(0,y,label,fontsize=size,va='top')
positions=[[.12,.60,.325,.22],[.605,.60,.325,.22],[.12,.335,.325,.18],[.605,.335,.325,.18],[.31,.055,.40,.165]]
for i,(pos,fn) in enumerate(zip(positions,funcs)):
    ax=fig.add_axes(pos);fn(ax);actual=ax.get_position();fig.text(actual.x0-.065,actual.y1+.022,chr(66+i),fontsize=10,fontweight='bold')
save(fig,'Fig1_discovery_definition_reconstructed','composites')
print('FIG1_COMPLETE')
