"""Render revision figures from analyst-exported processed tables only.

This script does not fit models, calculate statistics, select thresholds, or
recompute intervals. Time / 365 is solely an axis-unit conversion.
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
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data/processed';OUT=ROOT/'figures/panels'
OUT.mkdir(parents=True,exist_ok=True)
MM=1/25.4;LOW='#0072B2';HIGH='#D55E00';GRAY='#333333';TEAL='#009E73';PURPLE='#CC79A7'
mpl.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':8,
 'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.75,
 'xtick.major.width':.75,'ytick.major.width':.75,'lines.linewidth':1,'pdf.fonttype':42,
 'ps.fonttype':42,'svg.fonttype':'none','figure.dpi':300,'savefig.dpi':300,
 'savefig.bbox':'tight','savefig.pad_inches':.02})
font=str(font_manager.findfont('Arial',fallback_to_default=True));records=[];sources={}
def read(name):
    p=DATA/name;sources[name]={'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
    return pd.read_csv(p)
def clean(ax,forest=False):
    ax.spines[['top','right']].set_visible(False);ax.tick_params(length=2,pad=2)
    if forest:ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
def export(fig,name,w,h,inputs):
    with mpl.rc_context({'savefig.bbox':None}):
        for ext in ['pdf','png']:fig.savefig(OUT/(name+'.'+ext),dpi=300,bbox_inches=None)
    records.append({'name':name,'width_mm':w,'height_mm':h,'inputs':inputs,'font':font,'min_font_pt':7})
    plt.close(fig)
def fmtp(p):return f'{p:.3f}' if p>=.001 else f'{p:.2e}'
def err(ax,row,y,color=GRAY,marker='o',field='estimate'):
    v=row[field];ax.errorbar(v,y,xerr=[[v-row.ci_lower],[row.ci_upper-v]],fmt=marker,
        color=color,markersize=4,capsize=2,elinewidth=.9,zorder=3)
MODELS=[('OSARS_frozen','OSARS',GRAY),('Hong2022_8gene_log','Hong 8-gene',TEAL),
        ('Ma2024_TR_OSRG_3gene_zscore','Ma 3-gene',PURPLE)]
metric=read('discrimination_metrics_with_CI.csv')

# Fig3 A: five independent outer-fold predictions and conditional mean CI.
fold=read('nested_outer_fold_Cindex_CI.csv');nested=read('nested_Cindex_summary_CI.csv')
fig,ax=plt.subplots(figsize=(85*MM,87*MM));fig.subplots_adjust(left=.34,right=.96,bottom=.16,top=.89)
for y,k in zip([6,5,4,3,2],range(1,6)):
    r=fold[(fold.predictor=='OOF')&(fold.outer_fold==k)].iloc[0];err(ax,r,y)
    ax.text(.955,y,f'{r.estimate:.3f}',ha='right',va='center')
r=nested[nested.predictor=='Nested CV: mean outer-fold C-index'].iloc[0];err(ax,r,.7,TEAL,'D')
ax.text(.955,.7,f'{r.estimate:.3f}',ha='right',va='center')
r=metric[(metric.cohort=='ICGC-LIRI')&(metric.predictor=='new_TCGA_only_model')&(metric.metric=='C-index')].iloc[0]
err(ax,r,-.7,TEAL,'s');ax.text(.955,-.7,f'{r.estimate:.3f}',ha='right',va='center')
ax.axhline(1.35,color='.8',lw=.6);ax.axhline(0,color='.8',lw=.6)
ax.set_yticks([6,5,4,3,2,.7,-.7],['Outer fold 1','Outer fold 2','Outer fold 3','Outer fold 4','Outer fold 5','Outer mean','ICGC (n = 243)'])
ax.set_xlim(.4,.97);ax.set_ylim(-1.4,6.7);ax.set_xticks([.4,.5,.6,.7,.8,.9])
ax.axvline(.5,color='.7',lw=.65,ls=(0,(3,3)));ax.set_xlabel('Harrell’s C-index')
fig.text(.56,.95,'New TCGA training strategy',ha='center',fontsize=8)
clean(ax,True);export(fig,'fig3e',85,87,['nested_outer_fold_Cindex_CI.csv','nested_Cindex_summary_CI.csv','discrimination_metrics_with_CI.csv'])

# Fig3 B: published scores on identical complete cases, no direction changes.
fig,ax=plt.subplots(figsize=(85*MM,87*MM));fig.subplots_adjust(left=.35,right=.96,bottom=.16,top=.89)
for co,ys in [('TCGA-LIHC',[6,5,4]),('ICGC-LIRI_benchmark_completecases',[2,1,0])]:
    for y,(m,label,c) in zip(ys,MODELS):
        r=metric[(metric.cohort==co)&(metric.predictor==m)&(metric.metric=='C-index')].iloc[0]
        err(ax,r,y,c);ax.text(.968,y,f'{r.estimate:.3f}',ha='right',va='center')
ax.set_yticks([6,5,4,2,1,0],[m[1] for m in MODELS]*2)
ax.set_xlim(.5,.98);ax.set_ylim(-.7,7.2);ax.set_xticks([.5,.6,.7,.8,.9])
ax.set_xlabel('Harrell’s C-index');ax.axvline(.5,color='.7',lw=.65,ls=(0,(3,3)))
ax.text(.01,6.85,'TCGA-LIHC: apparent (n = 343)',transform=blended_transform_factory(ax.transAxes,ax.transData),fontsize=8)
ax.text(.01,2.85,'ICGC: retrospective (n = 141)',transform=blended_transform_factory(ax.transAxes,ax.transData),fontsize=8)
clean(ax,True);export(fig,'fig3f',85,87,['discrimination_metrics_with_CI.csv'])

# Fig3 C,D: 1- and 3-year cumulative/dynamic IPCW AUC.
for letter,co,title in [('c','TCGA-LIHC','TCGA-LIHC (n = 343)'),('d','ICGC-LIRI_benchmark_completecases','ICGC-LIRI (n = 141)')]:
    fig,ax=plt.subplots(figsize=(85*MM,70*MM));fig.subplots_adjust(left=.2,right=.96,bottom=.20,top=.74)
    for j,(m,label,c) in enumerate(MODELS):
        d=metric[(metric.cohort==co)&(metric.predictor==m)&(metric.metric=='AUC')&(metric.time_days.isin([365,1095]))].sort_values('time_days')
        x=np.array([1,3])+(j-1)*.10
        ax.errorbar(x,d.estimate,yerr=[d.estimate-d.ci_lower,d.ci_upper-d.estimate],fmt='o-',color=c,capsize=2,markersize=3,label=label)
    ax.set_xlim(.65,3.35);ax.set_ylim(.4,1);ax.set_xticks([1,3]);ax.set_yticks([.4,.6,.8,1])
    ax.axhline(.5,color='.7',lw=.65,ls=(0,(3,3)));ax.set_xlabel('Time (years)');ax.set_ylabel('Time-dependent AUC')
    fig.text(.57,.93,title,ha='center',fontsize=8)
    ax.legend(loc='lower center',bbox_to_anchor=(.5,1.05),ncol=3,frameon=False,columnspacing=.8,handlelength=1.0,handletextpad=.4)
    clean(ax);export(fig,'sfig_auc_'+{'c':'a','d':'b'}[letter],85,70,['discrimination_metrics_with_CI.csv'])

# Supplement: program direction and partial rank associations.
corr=read('OSARS_program_correlations.csv');partial=read('MP4_partial_rank_correlations.csv')
programs=[('ROS_mean_rank','ROS',GRAY),('MP4_mean_rank','MP4',TEAL),('MP4_without_ROS_mean_rank','MP4 minus ROS',PURPLE)]
for letter,co in [('a','TCGA-LIHC'),('b','ICGC-LIRI')]:
    d=corr[corr.cohort==co];fig,ax=plt.subplots(figsize=(85*MM,62*MM));fig.subplots_adjust(left=.37,right=.96,bottom=.24,top=.79)
    for y,(m,label,c) in zip([2,1,0],programs):
        r=d[d.program==m].iloc[0];err(ax,r,y,c,field='spearman_rho')
        ax.text(r.spearman_rho,y+.26,f'{r.spearman_rho:.3f}',ha='center',fontsize=7)
    ax.set_yticks([2,1,0],[x[1] for x in programs]);ax.set_xlim(-.75,.45);ax.set_xticks([-.6,-.3,0,.3]);ax.set_ylim(-.6,2.7)
    ax.axvline(0,color='.7',ls=(0,(3,3)),lw=.65);ax.set_xlabel('Spearman ρ with OSARS')
    fig.text(.6,.91,f'{co} (n = {int(d.iloc[0].n)})',ha='center',fontsize=8)
    clean(ax,True);export(fig,'fig3'+letter,85,62,['OSARS_program_correlations.csv'])
for letter,co in [('c','TCGA-LIHC'),('d','ICGC-LIRI')]:
    d=partial[partial.cohort==co];adjustments=list(d.adjustment.unique());h=80
    fig,ax=plt.subplots(figsize=(85*MM,h*MM));fig.subplots_adjust(left=.39,right=.96,bottom=.28,top=.87)
    for j,m in enumerate(['MP4_mean_rank','MP4_without_ROS_mean_rank']):
        for y,a in zip(range(len(adjustments)-1,-1,-1),adjustments):
            r=d[(d.program==m)&(d.adjustment==a)].iloc[0];err(ax,r,y+(j-.5)*.20,[TEAL,PURPLE][j],field='partial_spearman_rho')
    labs=[{'unadjusted':'Unadjusted','ROS':'ROS','ROS + stage':'ROS + stage','ROS + stage + ESTIMATE purity':'ROS + stage\n+ purity'}[a] for a in adjustments]
    ax.set_yticks(range(len(adjustments)-1,-1,-1),labs);ax.set_xlim(-.75,0);ax.set_xticks([-.6,-.4,-.2,0]);ax.set_ylim(-.6,len(adjustments)-.35)
    ax.axvline(0,color='.7',ls=(0,(3,3)),lw=.65);ax.set_xlabel('Partial rank correlation with OSARS')
    fig.text(.6,.94,co,ha='center',fontsize=8)
    handles=[Line2D([0],[0],color=TEAL,marker='o',lw=1,label='MP4'),Line2D([0],[0],color=PURPLE,marker='o',lw=1,label='MP4 minus ROS')]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.56,.02),ncol=2,frameon=False,columnspacing=1,handlelength=1)
    clean(ax,True);export(fig,'fig3'+letter,85,h,['MP4_partial_rank_correlations.csv'])

# Fig4 A-C and supporting KM curves; all risk tables/CI/censor coordinates supplied.
km=read('KM_coordinates.csv');risk=read('KM_risk_tables.csv');ktest=read('KM_logrank_tests.csv');cutoff=read('risk_cutoff_manifest.csv')
pkm=read('PDC_verified158_KM_coordinates.csv');prisk=read('PDC_verified158_KM_risk_tables.csv');ptest=read('PDC_verified158_KM_logrank_tests.csv');pcut=read('PDC_verified158_cutoff_manifest.csv')
def km_plot(co,predictor,title,name,w=170/3,h=95,endpoint='Overall survival',horizon=3,grouping='High'):
    is_pdc=co.startswith('PDC'); k=pkm if is_pdc else km;rr=prisk if is_pdc else risk;tt=ptest if is_pdc else ktest;cc=pcut if is_pdc else cutoff
    stat=tt[(tt.cohort==co)&(tt.grouping==grouping)].iloc[0];cut=cc[(cc.cohort==co)&(cc.predictor==predictor)].iloc[0]
    fig=plt.figure(figsize=(w*MM,h*MM));ax=fig.add_axes([.23,.39,.70,.46]);ra=fig.add_axes([.23,.06,.70,.19])
    for g,c in [('Low',LOW),('High',HIGH)]:
        d=k[(k.cohort==co)&(k.predictor==predictor)&(k.group==g)].sort_values('time_days')
        ax.step(d.time_days/365,d.survival,where='post',color=c,lw=1,label=g)
        ax.fill_between(d.time_days/365,d.ci_lower,d.ci_upper,step='post',color=c,alpha=.14,linewidth=0)
        v=d[d.n_censor>0];ax.scatter(v.time_days/365,v.survival,marker='|',s=9,color=c,linewidths=.6,zorder=3)
        r=rr[(rr.cohort==co)&(rr.predictor==predictor)&(rr.group==g)&(rr.time_days<=horizon*365)]
        y=1 if g=='Low' else 0
        for _,row in r.iterrows():ra.text(row.time_days/365,y,str(int(row.n_risk)),ha='center',va='center',color=c,fontsize=7)
        ra.text(-.10,y,g,transform=blended_transform_factory(ra.transAxes,ra.transData),ha='right',va='center',color=c,fontsize=7)
    ax.set_xlim(0,horizon);ax.set_ylim(0,1.025);ax.set_xticks(range(horizon+1));ax.set_yticks([0,.25,.5,.75,1]);ax.set_xlabel('Time (years)',labelpad=3);ax.set_ylabel(endpoint,labelpad=3)
    ax.legend(loc='lower left',frameon=False,handlelength=1.0,borderaxespad=.2,labelspacing=.3,bbox_to_anchor=(0,.12))
    ax.text(.02,.035,f'Log-rank P = {fmtp(stat.p)}',transform=ax.transAxes,fontsize=7)
    fig.text(.58,.98 if '\n' in title else .945,title,ha='center',va='top' if '\n' in title else 'baseline',fontsize=8)
    fig.text(.58,.887 if '\n' in title else .902,f'n = {int(stat.n)}; events = {int(stat.events)}',ha='center',fontsize=7)
    clean(ax);ra.set_xlim(0,horizon);ra.set_ylim(-.5,1.6);ra.axis('off');fig.text(.23,.287,'Number at risk',fontsize=7)
    inputs=['PDC_verified158_KM_coordinates.csv','PDC_verified158_KM_risk_tables.csv','PDC_verified158_KM_logrank_tests.csv','PDC_verified158_cutoff_manifest.csv'] if is_pdc else ['KM_coordinates.csv','KM_risk_tables.csv','KM_logrank_tests.csv','risk_cutoff_manifest.csv']
    export(fig,name,w,h,inputs)
km_plot('TCGA-LIHC','OSARS_frozen','TCGA-LIHC (training)','fig4a')
km_plot('ICGC-LIRI','OSARS_frozen','ICGC-LIRI','fig4b')
km_plot('PDC000198_OS_verified158','OSARS_protein_surrogate','PDC000198\nprotein surrogate','fig4c')
km_plot('ICGC-LIRI','OSARS_frozen_training_cutoff','ICGC: TCGA cutoff','sfig_clinical_a',85,95,grouping='High_train_cutoff')
km_plot('PDC000198_RFS_verified158','OSARS_protein_surrogate','PDC000198 protein surrogate','sfig_clinical_b',85,95,endpoint='Recurrence-free survival')

# Fig4 D average adjusted OSARS associations; primary PDC is the verified 158 only.
cox=read('Cox_forest_long.csv');sens=read('Cox_sensitivity_forest_long.csv')
fig,ax=plt.subplots(figsize=(85*MM,68*MM));fig.subplots_adjust(left=.40,right=.97,bottom=.25,top=.82)
labels=[]
for y,co,d,label in [(2,'TCGA-LIHC',cox,'TCGA-LIHC\n(n = 321)'),(1,'ICGC-LIRI',cox,'ICGC-LIRI\n(n = 243)'),(0,'PDC000198_OS_verified158',sens,'Protein surrogate\nPDC (n = 158)')]:
    r=d[(d.cohort==co)&(d.model=='multivariable_primary')&(d.term=='OSARS_SD')].iloc[0]
    err(ax,r,y,field='HR');ax.text(r.HR,y+.31,f'{r.HR:.2f} [{r.ci_lower:.2f}, {r.ci_upper:.2f}]',ha='center',fontsize=7);labels.append(label)
ax.set_xscale('log');ax.set_xlim(.7,7);ax.set_xticks([1,2,4],['1','2','4']);ax.minorticks_off();ax.set_ylim(-.6,2.8);ax.set_yticks([2,1,0],labels)
ax.axvline(1,color='.7',ls=(0,(3,3)),lw=.65);ax.set_xlabel('Adjusted hazard ratio per 1 SD')
fig.text(.60,.94,'Overall survival',ha='center',fontsize=8);clean(ax,True)
export(fig,'fig4d',85,68,['Cox_forest_long.csv','Cox_sensitivity_forest_long.csv'])

# Fig4 E clinical models trained on TCGA and evaluated in the existing ICGC cohort.
inc=read('clinical_incremental_discrimination_CI.csv');delta=read('clinical_incremental_paired_differences.csv')
fig,ax=plt.subplots(figsize=(85*MM,68*MM));fig.subplots_adjust(left=.37,right=.96,bottom=.31,top=.84)
models=[('Stage','Stage',GRAY),('Stage_plus_ROS','Stage + ROS',TEAL),('Stage_plus_OSARS','Stage + OSARS',PURPLE)]
for y,(m,label,c) in zip([2,1,0],models):
    r=inc[(inc.cohort=='ICGC-LIRI')&(inc.predictor==m)&(inc.metric=='C-index')].iloc[0]
    err(ax,r,y,c);ax.text(.96,y,f'{r.estimate:.3f}',ha='right',va='center')
di=delta[(delta.cohort=='ICGC-LIRI')&(delta.predictor=='Stage_plus_OSARS')&(delta.metric=='C-index')].iloc[0]
ax.set_yticks([2,1,0],[x[1] for x in models]);ax.set_ylim(-.6,2.7);ax.set_xlim(.5,.98);ax.set_xticks([.5,.6,.7,.8,.9]);ax.set_xlabel('Harrell’s C-index')
ax.axvline(.5,color='.7',ls=(0,(3,3)),lw=.65);fig.text(.59,.94,'ICGC-LIRI (n = 243)',ha='center',fontsize=8)
fig.text(.54,.06,f'ΔC = {di.delta:.3f} (95% CI {di.ci_lower:.3f} to {di.ci_upper:.3f})',ha='center',fontsize=7)
clean(ax,True);export(fig,'fig4e',85,68,['clinical_incremental_discrimination_CI.csv','clinical_incremental_paired_differences.csv'])

# Supplement clinical C: stage-stratified OS sensitivity and primary PDC RFS.
fig,ax=plt.subplots(figsize=(85*MM,65*MM));fig.subplots_adjust(left=.43,right=.97,bottom=.25,top=.82)
items=[('TCGA-LIHC','stage_stratified_sensitivity','TCGA\nstage-stratified OS'),('ICGC-LIRI','stage_stratified_sensitivity','ICGC\nstage-stratified OS'),('PDC000198_RFS_verified158','multivariable_primary','PDC protein\nadjusted RFS')]
for y,(co,mod,lab) in zip([2,1,0],items):
    r=sens[(sens.cohort==co)&(sens.model==mod)&(sens.term=='OSARS_SD')].iloc[0];err(ax,r,y,field='HR')
    ax.text(r.HR,y+.30,f'{r.HR:.2f} [{r.ci_lower:.2f}, {r.ci_upper:.2f}]',ha='center',fontsize=7)
ax.set_yticks([2,1,0],[x[2] for x in items]);ax.set_xscale('log');ax.set_xlim(.6,7);ax.set_xticks([1,2,4],['1','2','4']);ax.minorticks_off();ax.set_ylim(-.6,2.8)
ax.axvline(1,color='.7',ls=(0,(3,3)),lw=.65);ax.set_xlabel('Adjusted hazard ratio per 1 SD');clean(ax,True)
export(fig,'sfig_clinical_c',85,65,['Cox_sensitivity_forest_long.csv'])

# Supplement clinical D: TCGA time-interaction diagnostic, precomputed HR only.
tv=read('TCGA_OSARS_timevarying_HR.csv');fig,ax=plt.subplots(figsize=(85*MM,65*MM));fig.subplots_adjust(left=.2,right=.95,bottom=.23,top=.80)
ax.errorbar(tv.time_days/365,tv.HR_per_SD,yerr=[tv.HR_per_SD-tv.ci_lower,tv.ci_upper-tv.HR_per_SD],fmt='o-',color=GRAY,capsize=2,markersize=3)
ax.set_xlim(.25,3.25);ax.set_ylim(0,9);ax.set_xticks([.5,1,2,3]);ax.set_yticks([0,2,4,6,8]);ax.axhline(1,color='.7',lw=.65,ls=(0,(3,3)))
ax.set_xlabel('Time (years)');ax.set_ylabel('Adjusted HR per 1 SD');fig.text(.57,.93,'TCGA-LIHC: time-varying effect',ha='center',fontsize=8)
ax.text(.03,.94,f'Interaction P = {fmtp(tv.iloc[0].interaction_p)}',transform=ax.transAxes,va='top',fontsize=7)
clean(ax);export(fig,'sfig_clinical_d',85,65,['TCGA_OSARS_timevarying_HR.csv'])

# Supplement spline diagnostics: no re-fit, sampled curve and pointwise CI only.
spline=read('Cox_adjusted_spline_coordinates.csv');linearity=read('Cox_adjusted_linearity_diagnostics.csv')
for letter,co,title in [('a','TCGA-LIHC','TCGA-LIHC'),('b','ICGC-LIRI','ICGC-LIRI'),('c','PDC000198_OS_verified158','PDC000198 protein: OS')]:
    d=spline[spline.cohort==co];p=linearity[linearity.cohort==co].iloc[0]
    fig,ax=plt.subplots(figsize=(170/3*MM,80*MM));fig.subplots_adjust(left=.25,right=.94,bottom=.22,top=.83)
    ax.plot(d.OSARS_SD,d.HR,color=GRAY);ax.fill_between(d.OSARS_SD,d.ci_lower,d.ci_upper,color=GRAY,alpha=.13,linewidth=0)
    ax.set_yscale('log');ax.axhline(1,color='.7',lw=.65,ls=(0,(3,3)));ax.axvline(d.iloc[0].reference_OSARS_SD,color='.8',lw=.6,ls=(0,(2,3)))
    ticks={'TCGA-LIHC':[.01,.1,1,10,100],'ICGC-LIRI':[.2,.5,1,2,5,10],'PDC000198_OS_verified158':[.2,.5,1,2,5]}[co]
    ax.set_yticks(ticks,[str(t) for t in ticks]);ax.minorticks_off()
    ax.set_xlabel('OSARS (SD units)');ax.set_ylabel('Adjusted hazard ratio')
    fig.text(.59,.945,title,ha='center',fontsize=8);fig.text(.59,.885,f'Nonlinearity P = {fmtp(p.p)}',ha='center',fontsize=7)
    clean(ax);export(fig,'sfig_spline_'+letter,170/3,80,['Cox_adjusted_spline_coordinates.csv','Cox_adjusted_linearity_diagnostics.csv'])

(ROOT/'figures/panel_manifest.json').write_text(json.dumps({'panels':records,'sources':sources},indent=2))
print(f'Exported {len(records)} model validation panels.')
