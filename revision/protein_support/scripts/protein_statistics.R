a<-commandArgs(trailingOnly=TRUE);out<-a[1]
suppressPackageStartupMessages({library(survival);library(jsonlite)})
write_table<-function(x,nm)write.csv(x,file.path(out,paste0(nm,".csv")),row.names=FALSE,na="")
pairs<-read.csv(file.path(out,"nqo1_verified_pairs.csv"))
all_d<-read.csv(file.path(out,"nqo1_survival_input.csv"))
pairs$both_originally_measured<-tolower(as.character(pairs$both_originally_measured))=="true"
all_d$tissue_label_consistent<-tolower(as.character(all_d$tissue_label_consistent))=="true"
set.seed(20260905);B<-2000
pair_results<-list()
for(subset in c("All verified pairs","Both originally measured")) {
  pp<-if(subset=="All verified pairs")pairs else pairs[pairs$both_originally_measured,]
  dif<-pp$difference_tumor_minus_adjacent
  wt<-wilcox.test(pp$tumor,pp$adjacent,paired=TRUE,exact=FALSE,conf.int=TRUE)
  tt<-t.test(dif)
  boots<-replicate(B,{v<-sample(dif,replace=TRUE);c(mean(v),median(v),mean(v)/sd(v))})
  cis<-apply(boots,1,quantile,c(.025,.975),names=FALSE)
  pair_results[[subset]]<-data.frame(subset=subset,n_pairs=nrow(pp),
    tumor_median=median(pp$tumor),adjacent_median=median(pp$adjacent),
    mean_difference=mean(dif),mean_difference_lower=cis[1,1],mean_difference_upper=cis[2,1],
    median_difference=median(dif),median_difference_lower=cis[1,2],median_difference_upper=cis[2,2],
    paired_standardized_effect_dz=mean(dif)/sd(dif),dz_lower=cis[1,3],dz_upper=cis[2,3],
    hodges_lehmann_shift=unname(wt$estimate),hl_lower=wt$conf.int[1],hl_upper=wt$conf.int[2],
    wilcoxon_V=unname(wt$statistic),wilcoxon_p=wt$p.value,t_p=tt$p.value,
    positive_differences=sum(dif>0),negative_differences=sum(dif<0))
}
write_table(do.call(rbind,pair_results),"paired_NQO1_statistics")
cox_results<-list();ph<-list();linearity<-list();prediction<-list()
for(cohort in c("Tissue-label verified primary","Original source sensitivity")) {
d<-if(cohort=="Tissue-label verified primary")all_d[all_d$tissue_label_consistent,] else all_d
d$NQO1_z<-as.numeric(scale(d$NQO1))
d$male<-as.integer(tolower(d$Gender)=="male")
for(ep in c("OS","RFS")) {
  d$time<-if(ep=="OS")d$ostime else d$rsftime
  d$event<-if(ep=="OS")d$osevent else d$rsfevent
  stopifnot(all(d$time>0),all(d$event%in%c(0,1)))
  for(adjustment in c("Unadjusted","Age and sex adjusted")) {
    form<-if(adjustment=="Unadjusted")Surv(time,event)~NQO1_z else Surv(time,event)~NQO1_z+Age+male
    fit<-coxph(form,data=d,x=TRUE,y=TRUE);ss<-summary(fit)
    cox_results[[paste(cohort,ep,adjustment)]]<-data.frame(cohort=cohort,endpoint=ep,adjustment=adjustment,n=fit$n,events=fit$nevent,
      term=rownames(ss$coefficients),hr=ss$conf.int[,"exp(coef)"],lower=ss$conf.int[,"lower .95"],
      upper=ss$conf.int[,"upper .95"],p=ss$coefficients[,"Pr(>|z|)"],NQO1_unit="per 1 cohort SD")
    zz<-cox.zph(fit)$table
    ph[[paste(cohort,ep,adjustment)]]<-data.frame(cohort=cohort,endpoint=ep,adjustment=adjustment,term=rownames(zz),zz)
    flexible<-if(adjustment=="Unadjusted")coxph(Surv(time,event)~splines::ns(NQO1_z,df=3),data=d) else
      coxph(Surv(time,event)~splines::ns(NQO1_z,df=3)+Age+male,data=d)
    ll<-anova(fit,flexible,test="LRT")
    linearity[[paste(cohort,ep,adjustment)]]<-data.frame(cohort=cohort,endpoint=ep,adjustment=adjustment,nonlinearity_p=ll[2,"Pr(>|Chi|)"])
  }
}
}
cx<-do.call(rbind,cox_results)
cx$q_BH_primary_two_endpoints<-NA_real_
sel<-cx$term=="NQO1_z"&cx$adjustment=="Unadjusted"&cx$cohort=="Tissue-label verified primary"
cx$q_BH_primary_two_endpoints[sel]<-p.adjust(cx$p[sel],"BH")
write_table(cx,"NQO1_continuous_cox")
write_table(do.call(rbind,ph),"NQO1_cox_PH_diagnostics")
write_table(do.call(rbind,linearity),"NQO1_cox_linearity_diagnostics")
cor_results<-list()
d<-all_d[all_d$tissue_label_consistent,]
for(method in c("spearman","pearson")) {
  ct<-cor.test(d$NQO1,d$protein_OSARS_without_NQO1,method=method,exact=FALSE)
  boot<-replicate(B,{ii<-sample.int(nrow(d),nrow(d),replace=TRUE);cor(d$NQO1[ii],d$protein_OSARS_without_NQO1[ii],method=method)})
  ci<-quantile(boot,c(.025,.975),names=FALSE)
  cor_results[[method]]<-data.frame(method=method,n=nrow(d),estimate=unname(ct$estimate),
    lower=ci[1],upper=ci[2],p=ct$p.value,score="Existing protein OSARS surrogate with NQO1 contribution removed")
}
write_table(do.call(rbind,cor_results),"NQO1_surrogate_without_self_correlation")
write_table(d[,c("sample_id","case_id","NQO1","NQO1_z","protein_OSARS_without_NQO1")],"NQO1_surrogate_scatter_coordinates")
capture.output(sessionInfo(),file=file.path(out,"R_session_info.txt"))
cat("Verified pairs:",nrow(pairs),"; tumor clinical samples:",nrow(d),"; numerical exports complete.\n")
