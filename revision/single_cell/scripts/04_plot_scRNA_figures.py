"""Render frozen, processed single-cell results; no inferential calculations."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[1];P=ROOT/'data/processed'
OUT=ROOT/'figures';(OUT/'panels').mkdir(parents=True,exist_ok=True);(OUT/'composites').mkdir(exist_ok=True)
mpl.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':7,'axes.titlesize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.65,'xtick.major.width':.65,'ytick.major.width':.65,'lines.linewidth':.9,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.dpi':300})
LOW='#0072B2';HIGH='#D55E00';GRAY='#777777'
V=pd.read_csv(P/'validation_primary_scores.csv')
D=pd.read_csv(P/'discovery_validation_ready.csv')
ALL=pd.read_csv(P/'GSE149614_all_cells_embedding.csv',index_col=0)
M=pd.read_csv(P/'GSE149614_primary_malignant_embedding.csv',index_col=0)
M=M.join(V.set_index('cell')[['OS_composite','MP4_UCell','MP4_without_ROS_UCell']],validate='one_to_one')
assert len(M)==13691 and not M.OS_composite.isna().any()
T=pd.read_csv(P/'patient_level_program_tests.csv')
C=pd.read_csv(P/'within_patient_correlations.csv')
NMF=pd.read_csv(P/'cross_cohort_NMF_overlap.csv')
RNG=np.random.default_rng(20260905)

def tidy(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(length=2.4,pad=2)

def mapbase(ax):
    ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([])
    for s in ax.spines.values():s.set_visible(False)
    ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2')

def atlas(ax):
    palette={'Hepatocyte':HIGH,'T/NK':LOW,'Myeloid':'#E69F00','B':'#CC79A7','Endothelial':'#009E73','Fibroblast':'#777777'}
    for ct,col in palette.items():
        x=ALL.loc[ALL.celltype==ct]
        ax.scatter(x.UMAP1,x.UMAP2,s=.28,color=col,linewidths=0,rasterized=True)
    handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=col,markeredgecolor='none',markersize=3,label=ct) for ct,col in palette.items()]
    ax.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,-.07),ncol=3,frameon=False,columnspacing=.7,handletextpad=.3)
    mapbase(ax)

def feature(ax,field,label):
    x=M.sort_values(field)
    h=ax.scatter(x.UMAP1,x.UMAP2,c=x[field],s=.55,cmap='viridis',linewidths=0,rasterized=True)
    cb=ax.figure.colorbar(h,ax=ax,fraction=.045,pad=.025,shrink=.73)
    cb.set_label(label,fontsize=7);cb.ax.tick_params(labelsize=7,length=2)
    mapbase(ax)

def definitions(ax):
    defs=['OS_high_top20','OS_high_within_sample_top20','OS_UCell_high_within_sample_top20']
    labels=['Cohort-wide\ncomposite\nPrimary; n = 8','Within-patient\ncomposite; n = 10','Within-patient\nUCell; n = 10']
    for j,(prog,col,offset) in enumerate([('MP4',HIGH,.12),('MP4_without_ROS',LOW,-.12)]):
        q=T[(T.dataset=='GSE149614')&(T.program==prog)].set_index('definition').loc[defs]
        y=np.arange(3)+offset
        ax.errorbar(q.mean_difference,y,xerr=[q.mean_difference-q.ci_low,q.ci_high-q.mean_difference],fmt='o',color=col,ms=3,capsize=2,lw=.85)
    ax.axvline(0,color='#aaaaaa',ls='--',lw=.65)
    ax.set_yticks(range(3),labels);ax.invert_yaxis();ax.set_xlabel('Mean within-patient difference\n(OS-high − OS-low)')
    ax.legend(handles=[Line2D([],[],marker='o',color=HIGH,label='MP4',ms=3),Line2D([],[],marker='o',color=LOW,label='MP4 − ROS',ms=3)],loc='upper left',bbox_to_anchor=(0,1.17),ncol=2,frameon=False,handlelength=1,columnspacing=1)
    tidy(ax)

def correlations(ax):
    c=C[(C.dataset=='GSE149614')&(C.program=='MP4_without_ROS')].sort_values('sample')
    y=np.arange(len(c))
    for k,(_,r) in enumerate(c.iterrows()):ax.plot([r.rho,r.partial_rho_depth_mt],[k,k],color='#bbbbbb',lw=.8)
    ax.scatter(c.rho,y,s=10,c=HIGH,label='Unadjusted',zorder=3)
    ax.scatter(c.partial_rho_depth_mt,y,s=12,marker='s',facecolors='white',edgecolors=LOW,lw=.8,label='Depth/MT adjusted',zorder=3)
    ax.axvline(0,color='#aaaaaa',ls='--',lw=.65)
    ax.set_yticks(y,c['sample']);ax.invert_yaxis();ax.set_xlim(-1,1)
    ax.set_xticks([-1,-.5,0,.5,1]);ax.set_xlabel('Within-patient correlation (ρ)')
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.19),ncol=1,frameon=False,handlelength=1)
    tidy(ax)

def overlap(ax):
    n=NMF.pivot(index='discovery',columns='validation',values='jaccard').sort_index()
    h=ax.imshow(n,cmap='Blues',vmin=0,vmax=.35,aspect='equal')
    for i in range(4):
        for j in range(4):
            v=n.iloc[i,j]
            ax.text(j,i,f'{v:.2f}',ha='center',va='center',fontsize=7,color='white' if v>.19 else '#222222')
    ax.set_xticks(range(4),n.columns);ax.set_yticks(range(4),n.index)
    ax.set_xlabel('Validation NMF program');ax.set_ylabel('Discovery NMF program')
    for s in ax.spines.values():s.set_visible(False)
    cb=ax.figure.colorbar(h,ax=ax,fraction=.045,pad=.04,shrink=.7)
    cb.set_label('Jaccard index',fontsize=7);cb.ax.tick_params(labelsize=7,length=2)

def discovery_map(ax):
    x=D.sort_values('OS_composite')
    h=ax.scatter(x.UMAP1,x.UMAP2,c=x.OS_composite,cmap='viridis',s=.6,linewidths=0,rasterized=True)
    cb=ax.figure.colorbar(h,ax=ax,fraction=.045,pad=.025,shrink=.75);cb.set_label('OS composite')
    mapbase(ax)

def discovery_groups(ax):
    for flag,col,lab in [(False,LOW,'OS-low'),(True,HIGH,'OS-high')]:
        x=D[D.OS_high_top20==flag];ax.scatter(x.UMAP1,x.UMAP2,s=.6,color=col,linewidths=0,rasterized=True,label=lab)
    mapbase(ax);ax.legend(loc='upper center',bbox_to_anchor=(.5,-.06),ncol=2,frameon=False,markerscale=4)

def discovery_sample_fraction(ax):
    s=pd.read_csv(P/'discovery_sample_group_counts.csv')
    s=s.iloc[s['sample'].str.extract(r'(\d+)')[0].astype(int).argsort()]
    x=np.arange(len(s));ax.bar(x,s.fraction_high,color=HIGH,width=.7)
    ax.axhline(.2,color='#888888',ls='--',lw=.65)
    ax.set_xticks(x,s['sample'],rotation=45);ax.set_ylim(0,1)
    ax.set_ylabel('OS-high fraction');ax.set_xlabel('Tissue sample');tidy(ax)

def discovery_agreement(ax):
    s=pd.read_csv(P/'discovery_group_agreement.csv',index_col=0)
    h=ax.imshow(s,cmap='Blues',aspect='equal')
    for i in range(2):
        for j in range(2):ax.text(j,i,str(s.iloc[i,j]),ha='center',va='center',color='white' if s.iloc[i,j]>4000 else 'black')
    ax.set_xticks([0,1],['OS-low','OS-high']);ax.set_yticks([0,1],['OS-low','OS-high'])
    ax.set_xlabel('Revised top-20% definition');ax.set_ylabel('Historical saved definition')
    for sp in ax.spines.values():sp.set_visible(False)

def save(fig,name,folder):
    fig.savefig(OUT/folder/(name+'.pdf'),facecolor='white')
    fig.savefig(OUT/folder/(name+'.png'),dpi=300,facecolor='white')
    plt.close(fig)

funcs=[atlas,lambda ax:feature(ax,'OS_composite','OS composite'),lambda ax:feature(ax,'MP4_without_ROS_UCell','MP4 − ROS (UCell)'),definitions,correlations,overlap]
for i,fn in enumerate(funcs):
    fig,ax=plt.subplots(figsize=(3.35,2.9));fig.subplots_adjust(left=.38 if i==3 else .22,right=.94 if i==3 else .84,bottom=.25,top=.88);fn(ax)
    save(fig,'Fig2'+chr(65+i),'panels')
fig,axs=plt.subplots(3,2,figsize=(170/25.4,205/25.4))
fig.subplots_adjust(left=.12,right=.91,bottom=.09,top=.96,hspace=.63,wspace=.72)
for i,(ax,fn) in enumerate(zip(axs.flat,funcs)):
    fn(ax);pos=ax.get_position();fig.text(pos.x0-.062,pos.y1+.019,chr(65+i),fontsize=10,fontweight='bold')
save(fig,'Fig2_independent_single_cell_validation','composites')

funcs2=[discovery_map,discovery_groups,discovery_sample_fraction,discovery_agreement]
for i,fn in enumerate(funcs2):
    fig,ax=plt.subplots(figsize=(3.35,2.7));fig.subplots_adjust(left=.21,right=.85,bottom=.24,top=.91);fn(ax)
    save(fig,'FigS_discovery'+chr(65+i),'panels')
fig,axs=plt.subplots(2,2,figsize=(170/25.4,135/25.4));fig.subplots_adjust(left=.12,right=.92,bottom=.13,top=.95,hspace=.53,wspace=.63)
for i,(ax,fn) in enumerate(zip(axs.flat,funcs2)):
    fn(ax);pos=ax.get_position();fig.text(pos.x0-.06,pos.y1+.015,chr(65+i),fontsize=10,fontweight='bold')
save(fig,'Supplementary_discovery_definition_audit','composites')
print('SCRNA_FIGURES_COMPLETE')
