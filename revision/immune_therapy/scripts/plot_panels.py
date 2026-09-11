"""Render figure panels exclusively from data/processed. No statistical tests here."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import MaxNLocator

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/processed';OUT=ROOT/'figures/panels'
LOW='#0072B2';HIGH='#D55E00';COLORS={'Low':LOW,'High':HIGH}
mpl.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':8,
 'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.75,
 'xtick.major.width':.75,'ytick.major.width':.75,'lines.linewidth':1,'pdf.fonttype':42,
 'ps.fonttype':42,'svg.fonttype':'none','figure.dpi':300,'savefig.dpi':300,
 'savefig.bbox':'tight','savefig.pad_inches':.02})
FONT=font_manager.findfont('Arial',fallback_to_default=True)
MM=1/25.4
manifest=[]

def read(name):return pd.read_csv(DATA/name)
def clean(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(length=2,pad=2)
def export(fig,name,sources,size):
    OUT.mkdir(parents=True,exist_ok=True)
    with mpl.rc_context({'savefig.bbox':None}):
        for ext in ['pdf','png']:fig.savefig(OUT/(name+'.'+ext),bbox_inches=None,dpi=300)
    manifest.append({'file':name,'width_mm':size[0],'height_mm':size[1],'sources':sources,'font':'Arial','font_file':FONT,'base_font_pt':7})
    plt.close(fig)
def p_label(p):return f'P = {p:.4f}' if p>=.001 else f'P = {p:.2g}'
def q_label(q):return f'FDR = {q:.3f}' if q>=.001 else 'FDR < 0.001'

# A. Point estimates plus bootstrap confidence intervals, all 11 frozen features.
ist=read('immune_statistics.csv').sort_values('Standardized_difference',ascending=False)
fig,ax=plt.subplots(figsize=(85*MM,77*MM));fig.subplots_adjust(left=.495,right=.95,bottom=.20,top=.94)
labels=[]
for y,(_,r) in enumerate(ist.iterrows()):
    color=HIGH if r.Standardized_difference>0 else LOW
    ax.errorbar(r.Standardized_difference,y,xerr=[[r.Standardized_difference-r.CI_low],[r.CI_high-r.Standardized_difference]],
                fmt='o',color=color,markersize=3,capsize=2,elinewidth=.8,capthick=.7)
    if r.Significance!='ns':ax.text(r.CI_high+.06,y,r.Significance,va='center',ha='left',fontsize=7)
    method={'ssGSEA':'ssGSEA','CIBERSORT':'CIB','ESTIMATE':'EST'}[r.Method]
    labels.append(f'{r.Feature} ({method})')
ax.set_yticks(range(len(ist)),labels);ax.invert_yaxis();ax.set_ylim(len(ist)-.45,-.55)
ax.axvline(0,color='.7',lw=.65,ls=(0,(3,3)));ax.set_xlim(-.68,1.02);ax.set_xticks([-.5,0,.5,1])
ax.set_xlabel('Standardized difference\n(High − Low)',labelpad=4)
clean(ax);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
export(fig,'fig5a',['immune_statistics.csv'],(85,77))

def bxp(ax,feature,long,summary,stat,ylabel=None,show_n=False):
    boxes=[]
    for g in ['Low','High']:
        r=summary[(summary.Feature==feature)&(summary.Risk_Group==g)].iloc[0]
        boxes.append({k:float(r[k]) for k in ['q1','med','q3','whislo','whishi']})
    bp=ax.bxp(boxes,positions=[0,1],widths=.5,patch_artist=True,showfliers=False,
              medianprops={'color':'#202020','linewidth':.9},boxprops={'linewidth':.7},
              whiskerprops={'linewidth':.7},capprops={'linewidth':.7})
    for i,(g,patch) in enumerate(zip(['Low','High'],bp['boxes'])):
        patch.set_facecolor(COLORS[g]);patch.set_alpha(.35);patch.set_edgecolor(COLORS[g])
        d=long[(long.Feature==feature)&(long.Risk_Group==g)]
        ax.scatter(i+d.Jitter,d.Value,s=2.7,c=COLORS[g],alpha=.35,edgecolors='none',zorder=2)
    r=stat[stat.Feature==feature].iloc[0]
    extrema=summary[summary.Feature==feature];vmin=extrema.vmin.min();vmax=extrema.vmax.max();ran=vmax-vmin
    ax.set_ylim(vmin-.08*ran,vmax+.32*ran)
    by=vmax+.10*ran
    ax.plot([0,0,1,1],[by-.02*ran,by,by,by-.02*ran],color='#303030',lw=.65)
    ax.text(.5,by+.03*ran,str(r.Significance),ha='center',va='bottom',fontsize=7)
    ax.text(.5,1.055,feature,transform=ax.transAxes,ha='center',va='bottom',fontsize=8)
    ticks=['Low','High'] if not show_n else [f'Low\n(n = {int(r.Low_N)})',f'High\n(n = {int(r.High_N)})']
    ax.set_xticks([0,1],ticks);ax.set_xlim(-.5,1.5);ax.yaxis.set_major_locator(MaxNLocator(4))
    if ylabel:ax.set_ylabel(ylabel,labelpad=2)
    clean(ax)

# B. All six archived TIDE quantities, with individual observations.
tl=read('tide_patient_values.csv');tb=read('tide_box_summary.csv');ts=read('tide_statistics.csv')
fig,axes=plt.subplots(3,2,figsize=(85*MM,125*MM));fig.subplots_adjust(left=.145,right=.96,bottom=.065,top=.93,wspace=.52,hspace=.70)
for ax,f in zip(axes.flat,['TIDE','Exclusion','Dysfunction','MDSC','CAF','TAM M2']):
    bxp(ax,f,tl,tb,ts,ylabel='Score' if ax in axes[:,0] else None)
export(fig,'fig5b',['tide_patient_values.csv','tide_box_summary.csv','tide_statistics.csv'],(85,125))

# C. Group proportions with Wilson intervals; TIDE labels remain predicted.
rs=read('tide_response_statistics.csv')
fig,ax=plt.subplots(figsize=(85*MM,48*MM));fig.subplots_adjust(left=.20,right=.96,bottom=.27,top=.77)
for y,g in [(1,'Low'),(0,'High')]:
    r=rs[rs.Risk_Group==g].iloc[0];v=r.Proportion*100
    ax.errorbar(v,y,xerr=[[v-r.CI_low*100],[r.CI_high*100-v]],fmt='o',color=COLORS[g],markersize=4,capsize=2)
    ax.text(r.CI_high*100+3,y,f'{int(r.R)}/{int(r.N)}',ha='left',va='center',fontsize=7)
ax.text(.5,1.20,p_label(rs.P_value.iloc[0]),transform=ax.transAxes,ha='center',fontsize=7)
ax.set_yticks([1,0],['Low','High']);ax.set_ylim(-.65,1.65);ax.set_xlim(0,100);ax.set_xticks([0,25,50,75,100])
ax.set_xlabel('Predicted ICB response (%)',labelpad=4);ax.set_ylabel('OSARS',labelpad=4);clean(ax)
export(fig,'fig5c',['tide_response_statistics.csv'],(85,48))

# Treatment contexts are separate exploratory supplements; no effect is called efficacy.
tx=read('therapy_response_statistics.csv')
contexts=[('TACE (GSE104580)','TACE response (%)','TACE\nGSE104580','sfig_therapy_a'),
          ('Adjuvant sorafenib (GSE109211)','RFS-benefit signature class (%)','Adjuvant sorafenib\nGSE109211','sfig_therapy_b'),
          ('Melanoma TIL-ACT (GSE100797)','RECIST response (%)','Melanoma TIL-ACT\nGSE100797','sfig_therapy_c')]
for cohort,ylabel,cohort_label,name in contexts:
    d=tx[tx.Cohort==cohort]
    fig,ax=plt.subplots(figsize=(56.6666667*MM,68*MM));fig.subplots_adjust(left=.26,right=.95,bottom=.16,top=.77)
    for x,g in enumerate(['Low','High']):
        r=d[d.Risk_Group==g].iloc[0];v=r.Proportion*100
        ax.errorbar(x,v,yerr=[[v-r.CI_low*100],[r.CI_high*100-v]],fmt='o',color=COLORS[g],markersize=4,capsize=2)
        ax.text(x,r.CI_high*100+4,f'{int(r.R)}/{int(r.N)}',ha='center',va='bottom',fontsize=7)
    ax.set_xticks([0,1],['Low','High']);ax.set_xlim(-.6,1.6);ax.set_ylim(0,112);ax.set_yticks([0,25,50,75,100])
    ax.text(.5,1.035,p_label(d.P_value.iloc[0]),transform=ax.transAxes,ha='center',va='bottom',fontsize=7)
    ax.text(.5,1.19,cohort_label,transform=ax.transAxes,ha='center',va='bottom',fontsize=8)
    ax.set_ylabel(ylabel,labelpad=2);ax.set_xlabel('OSARS',labelpad=3);clean(ax)
    export(fig,name,['therapy_response_statistics.csv'],(56.6666667,68))

# DrugReflector: top 20 raw ranks are all retained, including unnamed compounds.
drug=read('drugreflector_top20.csv').sort_values('Rank_1based')
fig,ax=plt.subplots(figsize=(85*MM,110*MM));fig.subplots_adjust(left=.40,right=.96,bottom=.13,top=.96)
y=np.arange(len(drug));ax.hlines(y,0,drug.Model_probability_percent,color='#CACACA',lw=.7)
ax.scatter(drug.Model_probability_percent,y,color='#6B5D87',s=12,edgecolor='none',zorder=3)
ax.set_yticks(y,drug.Label);ax.invert_yaxis();ax.set_xlim(0,5);ax.set_ylim(len(drug)-.35,-.65)
ax.set_xlabel('DrugReflector model\nprobability (%)',labelpad=5)
clean(ax);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
export(fig,'sfig_drug_a',['drugreflector_top20.csv'],(85,110))

ck=read('checkpoint_patient_values.csv');cb=read('checkpoint_box_summary.csv');cs=read('checkpoint_statistics.csv')
fig,axes=plt.subplots(2,4,figsize=(170*MM,95*MM));fig.subplots_adjust(left=.065,right=.98,bottom=.08,top=.90,wspace=.65,hspace=.7)
for ax,g in zip(axes.flat,['PDCD1','CD274','CTLA4','LAG3','TIGIT','HAVCR2','VSIR','IDO1']):
    bxp(ax,g,ck,cb,cs,ylabel='Expression (source units)' if ax in axes[:,0] else None)
export(fig,'sfig_checkpoints_a',['checkpoint_patient_values.csv','checkpoint_box_summary.csv','checkpoint_statistics.csv'],(170,95))

# Direct fit-quality sensitivity for the two strongest archived CIBERSORT signals.
cib=read('cibersort_22_fraction_sensitivity.csv')
for feature,label,stem in [('T cells regulatory (Tregs)','Tregs','sfig_cibersort_a'),('Macrophages M0','Macrophages M0','sfig_cibersort_b')]:
    d=cib[cib.Feature==feature]
    fig,ax=plt.subplots(figsize=(85*MM,63*MM));fig.subplots_adjust(left=.43,right=.96,bottom=.24,top=.80)
    labels=[]
    for y,flt in [(1,'all_343'),(0,'CIBERSORT_P_lt_0.05')]:
        r=d[d.Filter==flt].iloc[0];v=r.Difference*100
        ax.errorbar(v,y,xerr=[[v-r.Difference_CI_low*100],[r.Difference_CI_high*100-v]],fmt='o',color='#333333',
                    markerfacecolor='#333333' if y else 'white',markersize=4,capsize=2,lw=.9)
        ax.text(.98,(y+.45)/2.25,q_label(r.FDR),transform=ax.transAxes,ha='right',va='bottom',fontsize=7)
        name='All estimates' if y else 'CIBERSORT P < 0.05'
        labels.append((y,f'{name}\nLow {int(r.Low_N)} | High {int(r.High_N)}'))
    ax.set_yticks([x[0] for x in labels],[x[1] for x in labels]);ax.set_ylim(-.65,1.6)
    ax.axvline(0,color='.65',lw=.65,ls=(0,(3,3)))
    ax.set_xlim(-12,16);ax.set_xticks([-10,0,10])
    assert (d.Difference_CI_low*100>=-12).all() and (d.Difference_CI_high*100<=16).all()
    ax.set_xlabel('Difference (percentage points)\n(High − Low)',labelpad=4)
    ax.text(.5,1.13,label,transform=ax.transAxes,ha='center',va='bottom',fontsize=8)
    clean(ax);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
    export(fig,stem,['cibersort_22_fraction_sensitivity.csv'],(85,63))
(ROOT/'qc/panel_export_manifest.json').write_text(json.dumps(manifest,indent=2))
print(f'Exported {len(manifest)} panels as editable-font PDF and 300 dpi PNG.')
