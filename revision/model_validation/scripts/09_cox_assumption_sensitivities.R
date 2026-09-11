suppressPackageStartupMessages(library(survival))
script_arg<-grep('^--file=',commandArgs(trailingOnly=FALSE),value=TRUE)[1]
module_dir<-normalizePath(file.path(dirname(sub('^--file=','',script_arg)),'..'),mustWork=FALSE)
out<-file.path(module_dir,'data/processed')
cohort_names<-c('TCGA-LIHC','ICGC-LIRI','ICGC-LIRI_unique_profiles','PDC000198_OS_verified158','PDC000198_RFS_verified158')
forest<-list();ph<-list();spline_curves<-list();lr<-list();time_curves<-list()
export_fit<-function(fit,co,name){
 s<-summary(fit);b<-s$coefficients;ci<-s$conf.int
 forest[[paste(co,name)]]<<-data.frame(cohort=co,model=name,term=rownames(b),beta=b[,1],SE=b[,3],HR=ci[,1],ci_lower=ci[,3],ci_upper=ci[,4],p=b[,5],n=fit$n,events=fit$nevent,formula=paste(deparse(fit$formula),collapse=' '))
 z<-tryCatch(cox.zph(fit)$table,error=function(e)NULL)
 if(!is.null(z))ph[[paste(co,name)]]<<-data.frame(cohort=co,model=name,term=rownames(z),chisq=z[,1],df=z[,2],p=z[,3])
}
for(co in cohort_names){
 protein<-grepl('PDC',co)
 d<-read.csv(file.path(out,paste0(co,if(protein)'_all_scores.csv'else'_clinical_analysis_rows.csv')))
 covars<-if(protein)c('age10','male')else c('age10','male','stage_advanced')
 d<-d[complete.cases(d[,c('time','event','OSARS_SD',covars)]),]
 form<-as.formula(paste('Surv(time,event)~OSARS_SD+',paste(covars,collapse='+')))
 flin<-coxph(form,data=d,x=TRUE,y=TRUE)
 if(protein)export_fit(flin,co,'multivariable_primary')
 fgrp<-coxph(as.formula(paste('Surv(time,event)~High+',paste(covars,collapse='+'))),data=d,x=TRUE,y=TRUE)
 export_fit(fgrp,co,'multivariable_High')
 if(!protein){
  # Stage violates PH in ICGC; allow separate baseline hazards by early/advanced stage.
  fs<-coxph(Surv(time,event)~OSARS_SD+age10+male+strata(stage_advanced),data=d,x=TRUE,y=TRUE)
  export_fit(fs,co,'stage_stratified_sensitivity')
 }
 # Covariate-adjusted natural spline, fixed 3 df, ref at complete-case median.
 basis<-splines::ns(d$OSARS_SD,df=3)
 d$b1<-basis[,1];d$b2<-basis[,2];d$b3<-basis[,3]
 fsp<-coxph(as.formula(paste('Surv(time,event)~b1+b2+b3+',paste(covars,collapse='+'))),data=d,x=TRUE,y=TRUE)
 stat<-2*(fsp$loglik[2]-flin$loglik[2])
 lr[[co]]<-data.frame(cohort=co,n=nrow(d),events=sum(d$event),LR_chisq=stat,df=2,p=pchisq(stat,2,lower.tail=FALSE),test='adjusted 3-df natural spline versus adjusted linear OSARS; exploratory diagnostic')
 xx<-seq(quantile(d$OSARS_SD,.01),quantile(d$OSARS_SD,.99),length.out=101);ref<-median(d$OSARS_SD)
 delta<-sweep(predict(basis,xx),2,as.numeric(predict(basis,ref)),'-')
 bb<-coef(fsp)[c('b1','b2','b3')];V<-vcov(fsp)[c('b1','b2','b3'),c('b1','b2','b3')]
 lp<-as.numeric(delta%*%bb);se<-sqrt(rowSums((delta%*%V)*delta))
 spline_curves[[co]]<-data.frame(cohort=co,OSARS_SD=xx,reference_OSARS_SD=ref,HR=exp(lp),ci_lower=exp(lp-1.96*se),ci_upper=exp(lp+1.96*se),n=nrow(d),events=sum(d$event),covariates=paste(covars,collapse=';'))
 if(co=='TCGA-LIHC'){
  # Address the score-specific PH diagnostic with a log-time interaction.
  ft<-coxph(Surv(time,event)~OSARS_SD+tt(OSARS_SD)+age10+male+stage_advanced,data=d,
    tt=function(x,t,...)x*log(pmax(t,1)/365),x=TRUE,y=TRUE)
  export_fit(ft,co,'OSARS_logtime_interaction_sensitivity')
  tt<-c(182.5,365,730,1095);M<-cbind(1,log(tt/365));b<-coef(ft)[c('OSARS_SD','tt(OSARS_SD)')];v<-vcov(ft)[c('OSARS_SD','tt(OSARS_SD)'),c('OSARS_SD','tt(OSARS_SD)')]
  lp<-as.numeric(M%*%b);se<-sqrt(rowSums((M%*%v)*M))
  time_curves[[co]]<-data.frame(cohort=co,time_days=tt,HR_per_SD=exp(lp),ci_lower=exp(lp-1.96*se),ci_upper=exp(lp+1.96*se),interaction_p=summary(ft)$coefficients['tt(OSARS_SD)','Pr(>|z|)'])
 }
}
write.csv(do.call(rbind,forest),file.path(out,'Cox_sensitivity_forest_long.csv'),row.names=FALSE)
write.csv(do.call(rbind,ph),file.path(out,'Cox_sensitivity_PH_diagnostics.csv'),row.names=FALSE)
write.csv(do.call(rbind,lr),file.path(out,'Cox_adjusted_linearity_diagnostics.csv'),row.names=FALSE)
write.csv(do.call(rbind,spline_curves),file.path(out,'Cox_adjusted_spline_coordinates.csv'),row.names=FALSE)
write.csv(do.call(rbind,time_curves),file.path(out,'TCGA_OSARS_timevarying_HR.csv'),row.names=FALSE)
