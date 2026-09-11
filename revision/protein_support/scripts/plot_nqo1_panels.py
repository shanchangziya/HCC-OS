"""Plot orthogonal NQO1 evidence from verified, analyst-processed tables only."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter

ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data/processed';OUT=ROOT/'figures/panels'
OUT.mkdir(parents=True,exist_ok=True);MM=1/25.4;W=170/3
mpl.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':8,
 'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.75,
 'xtick.major.width':.75,'ytick.major.width':.75,'lines.linewidth':1,'pdf.fonttype':42,
 'ps.fonttype':42,'svg.fonttype':'none','figure.dpi':300,'savefig.dpi':300,'savefig.bbox':'tight','savefig.pad_inches':.02})
font=str(font_manager.findfont('Arial',fallback_to_default=True));source=[];panels=[]
def read(name):
    p=DATA/name;source.append({'file':name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()});return pd.read_csv(p)
def clean(ax):ax.spines[['top','right']].set_visible(False);ax.tick_params(length=2,pad=2)
def save(fig,name,files):
    with mpl.rc_context({'savefig.bbox':None}):
        for ext in ['pdf','png']:fig.savefig(OUT/(name+'.'+ext),bbox_inches=None,dpi=300)
    panels.append({'name':name,'width_mm':W,'height_mm':85,'files':files});plt.close(fig)

pairs=read('nqo1_verified_pairs.csv');ps=read('paired_NQO1_statistics.csv')
s=ps[ps.subset=='All verified pairs'].iloc[0]
assert len(pairs)==s.n_pairs and pairs.case_id.is_unique and pairs.both_originally_measured.all()
fig,ax=plt.subplots(figsize=(W*MM,85*MM));fig.subplots_adjust(left=.26,right=.96,bottom=.18,top=.79)
for _,r in pairs.iterrows():ax.plot([0,1],[r.adjacent,r.tumor],color='#AAAAAA',lw=.4,alpha=.25,zorder=1)
ax.scatter(np.zeros(len(pairs)),pairs.adjacent,s=4,color='#0072B2',alpha=.65,linewidths=0,zorder=2)
ax.scatter(np.ones(len(pairs)),pairs.tumor,s=4,color='#D55E00',alpha=.65,linewidths=0,zorder=2)
ax.hlines([s.adjacent_median,s.tumor_median],[-.18,.82],[.18,1.18],color='#222222',lw=1.3,zorder=3)
ax.set_xlim(-.5,1.5);ax.set_xticks([0,1],['Adjacent','Tumor'])
ax.set_ylabel('NQO1 abundance\n(normalized log2 ratio)',labelpad=3)
ax.text(.5,1.16,f'n = {int(s.n_pairs)} pairs',transform=ax.transAxes,ha='center',fontsize=7)
ax.text(.5,1.07,f'P = {s.wilcoxon_p:.2e}',transform=ax.transAxes,ha='center',fontsize=7)
clean(ax);save(fig,'sfig_nqo1_a',['nqo1_verified_pairs.csv','paired_NQO1_statistics.csv'])

cox=read('NQO1_continuous_cox.csv')
c=cox[(cox.cohort=='Tissue-label verified primary')&(cox.term=='NQO1_z')]
assert len(c)==4 and (c.n==158).all()
fig,ax=plt.subplots(figsize=(W*MM,85*MM));fig.subplots_adjust(left=.43,right=.95,bottom=.18,top=.79)
labels=[]
for y,(_,r) in zip([3,2,1,0],c.iterrows()):
    adjusted=r.adjustment=='Age and sex adjusted'
    ax.errorbar(r.hr,y,xerr=[[r.hr-r.lower],[r.upper-r.hr]],fmt='o',color='#333333',
                markerfacecolor='white' if adjusted else '#333333',markersize=4,capsize=2,lw=.8)
    labels.append(f'{r.endpoint}\n'+('Age + sex' if adjusted else 'Unadjusted'))
ax.set_yticks([3,2,1,0],labels);ax.set_ylim(-.7,3.7);ax.set_xscale('log');ax.set_xlim(.65,1.5)
ax.xaxis.set_major_locator(FixedLocator([.75,1,1.25]));ax.xaxis.set_major_formatter(FixedFormatter(['0.75','1.00','1.25']));ax.xaxis.set_minor_formatter(NullFormatter());ax.minorticks_off()
ax.axvline(1,color='.65',lw=.65,ls=(0,(3,3)))
ax.set_xlabel('Hazard ratio\nper 1 SD NQO1',labelpad=4)
ax.text(.5,1.16,'n = 158',transform=ax.transAxes,ha='center',fontsize=7)
ax.text(.5,1.07,'OS / RFS',transform=ax.transAxes,ha='center',fontsize=7)
clean(ax);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
save(fig,'sfig_nqo1_b',['NQO1_continuous_cox.csv'])

sc=read('NQO1_surrogate_scatter_coordinates.csv');cs=read('NQO1_surrogate_without_self_correlation.csv')
r=cs[cs.method=='spearman'].iloc[0];assert len(sc)==r.n
fig,ax=plt.subplots(figsize=(W*MM,85*MM));fig.subplots_adjust(left=.29,right=.96,bottom=.18,top=.79)
ax.scatter(sc.NQO1,sc.protein_OSARS_without_NQO1,s=7,color='#6B5D87',alpha=.7,linewidths=0)
ax.set_xlabel('NQO1 abundance\n(normalized log2 ratio)',labelpad=4)
ax.set_ylabel('Protein OSARS surrogate\n(excluding NQO1)',labelpad=3)
ax.text(.5,1.16,f'Spearman ρ = {r.estimate:.3f}',transform=ax.transAxes,ha='center',fontsize=7)
ax.text(.5,1.07,f'P = {r.p:.4f}; n = {int(r.n)}',transform=ax.transAxes,ha='center',fontsize=7)
clean(ax);save(fig,'sfig_nqo1_c',['NQO1_surrogate_scatter_coordinates.csv','NQO1_surrogate_without_self_correlation.csv'])
(ROOT/'figures/panel_manifest.json').write_text(json.dumps({'panels':panels,'source':source,'font':font,'statistics_recomputed':False},indent=2))
print('Exported 3 NQO1 protein panels.')
