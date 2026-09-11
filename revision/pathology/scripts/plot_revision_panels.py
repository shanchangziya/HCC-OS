"""Render the corrected pathology supplement from analyst-exported processed data.
No patient filtering, statistical tests, model fitting, or CI calculation occurs here.
"""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.transforms import blended_transform_factory

ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data/processed';OUT=ROOT/'figures/panels'
OUT.mkdir(parents=True,exist_ok=True)
MM=1/25.4;COLORS={'Low':'#0072B2','High':'#D55E00'}
mpl.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':8,
 'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.75,
 'xtick.major.width':.75,'ytick.major.width':.75,'lines.linewidth':1,'pdf.fonttype':42,
 'ps.fonttype':42,'svg.fonttype':'none','figure.dpi':300,'savefig.dpi':300,
 'savefig.bbox':'tight','savefig.pad_inches':.02})
font=str(font_manager.findfont('Arial',fallback_to_default=True));records=[];sources=[]
def read(name):
    path=DATA/name;sources.append({'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    return pd.read_csv(path)
def clean(ax):
    ax.spines[['top','right']].set_visible(False);ax.tick_params(length=2,pad=2)
def export(fig,name,w,h,inputs):
    with mpl.rc_context({'savefig.bbox':None}):
        for ext in ['pdf','png']:fig.savefig(OUT/(name+'.'+ext),dpi=300,bbox_inches=None)
    records.append({'name':name,'width_mm':w,'height_mm':h,'inputs':inputs,'font':font,'min_font_pt':7})
    plt.close(fig)
def fmtp(p):
    if p>=.001:return f'{p:.4f}'
    return f'{p:.2e}'

km=read('km_coordinates.csv');censor=read('km_censor_coordinates.csv');risk=read('km_risk_table.csv')
groupstat=read('survival_group_statistics.csv');perf=read('performance_metrics.csv');delta=read('paired_incremental_performance.csv')
for split,name,label in [('Train','sfig_pathology_a','Training'),('Test','sfig_pathology_b','Historical test')]:
    stat=groupstat[groupstat.split==split].iloc[0]
    fig=plt.figure(figsize=(85*MM,95*MM))
    ax=fig.add_axes([.17,.39,.79,.50]);ra=fig.add_axes([.17,.075,.79,.185])
    for g in ['Low','High']:
        d=km[(km.split==split)&(km.group==g)].sort_values('time_years')
        ax.step(d.time_years,d.survival,where='post',color=COLORS[g],lw=1,label=f'{g} (n = {int(stat.n_low if g=="Low" else stat.n_high)})')
        ax.fill_between(d.time_years,d.lower,d.upper,step='post',color=COLORS[g],alpha=.14,linewidth=0)
        c=censor[(censor.split==split)&(censor.group==g)]
        ax.scatter(c.time_years,c.survival,marker='|',s=10,color=COLORS[g],linewidths=.65,zorder=3)
    ax.set_xlim(0,5);ax.set_ylim(0,1.025);ax.set_xticks(range(6));ax.set_yticks([0,.25,.5,.75,1])
    ax.set_xlabel('Time (years)',labelpad=4);ax.set_ylabel('Overall survival probability',labelpad=3)
    ax.text(.03,.13,f'HR = {stat.hr_high_vs_low:.2f} (95% CI {stat.hr_lower:.2f}–{stat.hr_upper:.2f})',transform=ax.transAxes,fontsize=7)
    ax.text(.03,.045,f'Log-rank P = {fmtp(stat.logrank_p)}',transform=ax.transAxes,fontsize=7)
    ax.legend(loc='upper right',frameon=False,handlelength=1.5,borderaxespad=.25,labelspacing=.35)
    ax.text(.5,1.14,f'{label} (n = {int(stat.n)})',transform=ax.transAxes,ha='center',fontsize=8)
    clean(ax)
    ra.set_xlim(0,5);ra.set_ylim(-.5,1.6);ra.axis('off')
    blend=blended_transform_factory(ra.transAxes,ra.transData)
    for y,g in [(1,'Low'),(0,'High')]:
        r=risk[(risk.split==split)&(risk.group==g)]
        for _,row in r.iterrows():ra.text(row.time_years,y,str(int(row.n_risk)),ha='center',va='center',fontsize=7,color=COLORS[g])
        ra.text(-.085,y,g,transform=blend,ha='right',va='center',color=COLORS[g],fontsize=7)
    fig.text(.17,.284,'Number at risk',fontsize=7)
    export(fig,name,85,95,['km_coordinates.csv','km_censor_coordinates.csv','km_risk_table.csv','survival_group_statistics.csv'])

# C. Unchanged frozen pathology predictor, apparent training and historical test.
p=perf[(perf.subset=='All')&(perf.model=='Pathology')&(perf.metric=='C_index')]
fig,ax=plt.subplots(figsize=(85*MM,60*MM));fig.subplots_adjust(left=.35,right=.96,bottom=.24,top=.81)
for y,split in [(1,'Train'),(0,'Test')]:
    r=p[p.split==split].iloc[0]
    ax.errorbar(r.estimate,y,xerr=[[r.estimate-r.lower],[r.upper-r.estimate]],fmt='o',color='#333333',markersize=4,capsize=2)
    ax.text(r.estimate,y+.22,f'{r.estimate:.3f} [{r.lower:.3f}, {r.upper:.3f}]',ha='center',va='bottom',fontsize=7,bbox={'facecolor':'white','edgecolor':'none','pad':.4})
ax.set_yticks([1,0],['Training\n(n = 230)','Historical test\n(n = 100)'])
ax.set_xlim(.4,.9);ax.set_ylim(-.55,1.65);ax.set_xticks([.4,.5,.6,.7,.8,.9])
ax.axvline(.5,color='.7',lw=.65,ls=(0,(3,3)))
ax.set_xlabel('Harrell’s C-index',labelpad=4)
ax.text(.5,1.14,'Pathology score',transform=ax.transAxes,ha='center',fontsize=8)
clean(ax);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
export(fig,'sfig_pathology_c',85,60,['performance_metrics.csv'])

# D. Paired clinical-complete test sample; all three predictors use the same 94 cases.
p=perf[(perf.split=='Test')&(perf.subset=='Clinical_complete')&(perf.metric=='C_index')]
di=delta[(delta.split=='Test')&(delta.metric=='C_index')].iloc[0]
fig,ax=plt.subplots(figsize=(85*MM,60*MM));fig.subplots_adjust(left=.37,right=.96,bottom=.27,top=.81)
models=[('Clinical','Clinical'),('Pathology','Pathology'),('Clinical_plus_pathology','Clinical +\npathology')]
for y,(model,label) in zip([2,1,0],models):
    r=p[p.model==model].iloc[0]
    ax.errorbar(r.estimate,y,xerr=[[r.estimate-r.lower],[r.upper-r.estimate]],fmt='o',color='#333333',markersize=4,capsize=2)
    ax.text(.89,y,f'{r.estimate:.3f}',ha='right',va='center',fontsize=7)
ax.set_yticks([2,1,0],[x[1] for x in models]);ax.set_xlim(.4,.9);ax.set_ylim(-.6,2.6);ax.set_xticks([.4,.5,.6,.7,.8,.9])
ax.axvline(.5,color='.7',lw=.65,ls=(0,(3,3)));ax.set_xlabel('Harrell’s C-index',labelpad=4)
ax.text(.5,1.14,'Historical test (n = 94)',transform=ax.transAxes,ha='center',fontsize=8)
fig.text(.52,.055,f'ΔC = {di.difference:.3f} (95% CI {di.lower:.3f} to {di.upper:.3f})',ha='center',fontsize=7)
clean(ax);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
export(fig,'sfig_pathology_d',85,60,['performance_metrics.csv','paired_incremental_performance.csv'])
(ROOT/'figures/panel_manifest.json').write_text(json.dumps({'panels':records,'sources':sources},indent=2))
print(f'Exported {len(records)} corrected pathology panels.')
