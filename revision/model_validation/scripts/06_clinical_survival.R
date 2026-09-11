suppressPackageStartupMessages(library(survival))
script_arg<-grep('^--file=',commandArgs(trailingOnly=FALSE),value=TRUE)[1]
out<-normalizePath(file.path(dirname(sub('^--file=','',script_arg)),'..'),mustWork=FALSE)
processed<-file.path(out,'data/processed');raw<-file.path(out,'data/raw')
study<-Sys.getenv('HCC_OS_LEGACY_PROJECT_ROOT',normalizePath(file.path(out,'../..'),mustWork=FALSE))
write_out<-function(x,n)write.csv(x,file.path(processed,n),row.names=FALSE)
risks<-function(co){x<-read.csv(file.path(processed,paste0(co,'_all_scores.csv')),check.names=FALSE);x}
tc<-risks('TCGA-LIHC');cl<-read.csv(file.path(raw,'TCGA_clinical_source.csv'))
extra_source<-Sys.getenv('HCC_OS_TCGA_CLINICAL_GRADE',file.path(study,'FIG4_clinical_utility/Fig4_TCGA_OSARS_clinical_all_matched.csv'))
extra<-read.csv(extra_source)
stopifnot(!anyDuplicated(tc$ID),!anyDuplicated(cl$ID),setequal(tc$ID,cl$ID))
tc<-merge(tc,cl[,c('ID','age','gender','stage')],by='ID',sort=FALSE)
tc<-merge(tc,extra[,c('ID','Grade')],by='ID',all.x=TRUE,sort=FALSE)
tc$stage_num<-match(tc$stage,c('I','II','III','IV'));tc$grade<-tc$Grade
ic<-risks('ICGC-LIRI');cl<-read.csv(file.path(raw,'ICGC_clinical_source.csv'))
stopifnot(setequal(ic$ID,cl$ID))
ic<-merge(ic,cl[,c('ID','age','gender','stage','grade')],by='ID',sort=FALSE)
ic$stage_num<-as.numeric(ic$stage)
cohorts<-list('TCGA-LIHC'=tc,'ICGC-LIRI'=ic)
dedup<-read.csv(file.path(processed,'ICGC_identical_expression_profile_audit.csv'))
cohorts[['ICGC-LIRI_unique_profiles']]<-ic[ic$ID %in% dedup$ID[tolower(as.character(dedup$retain_profile_sensitivity))=='true'],]
for(co in names(cohorts)){
 d<-cohorts[[co]];d$OSARS_SD<-(d$OSARS_frozen-mean(d$OSARS_frozen))/sd(d$OSARS_frozen)
 d$age10<-d$age/10;d$male<-ifelse(is.na(d$gender),NA,as.integer(tolower(d$gender)=='male'))
 d$stage_advanced<-ifelse(is.na(d$stage_num),NA,as.integer(d$stage_num>=3))
 d$grade_high<-ifelse(d$grade %in% c('G1','G2','I','II'),0,ifelse(d$grade %in% c('G3','G4','III','IV'),1,NA))
 d$High<-as.integer(d$OSARS_frozen>median(d$OSARS_frozen))
 d$High_train_cutoff<-as.integer(d$OSARS_frozen>median(tc$OSARS_frozen))
 cohorts[[co]]<-d
}
cl<-read.csv(file.path(raw,'CPTAC_clinical_source.csv'),check.names=FALSE)
for(endpoint in c('OS','RFS')){
 co<-paste0('CPTAC_',endpoint);d<-risks(co)
 j<-match(d$ID,cl$ID);stopifnot(!anyNA(j),all.equal(d$time,cl[[if(endpoint=='OS')'ostime'else'rsftime']][j]*365/12)==TRUE)
 d$age<-cl$Age[j];d$gender<-cl$Gender[j];d$age10<-d$age/10;d$male<-as.integer(tolower(d$gender)=='male')
 d$ALB<-cl$ALB[j];d$logAFP<-log2(cl$AFP[j]+1)
 d$OSARS_SD<-(d$OSARS_protein_surrogate-mean(d$OSARS_protein_surrogate))/sd(d$OSARS_protein_surrogate)
 d$High<-as.integer(d$OSARS_protein_surrogate>median(d$OSARS_protein_surrogate));cohorts[[co]]<-d
}
all_results<-list();ph<-list();descriptives<-list();nonlinear<-list();km_tests<-list()
fit_export<-function(d,co,name,terms){
 f<-as.formula(paste('Surv(time,event)~',paste(terms,collapse='+')))
 dd<-d[complete.cases(d[,c('time','event',terms)]),]
 fit<-coxph(f,data=dd,ties='efron',x=TRUE,y=TRUE,model=TRUE)
 s<-summary(fit);b<-s$coefficients;ci<-s$conf.int
 z<-data.frame(cohort=co,model=name,term=rownames(b),beta=b[,1],SE=b[,3],HR=ci[,1],ci_lower=ci[,3],ci_upper=ci[,4],p=b[,5],n=nrow(dd),events=sum(dd$event),formula=paste(deparse(f),collapse=' '),stringsAsFactors=FALSE)
 all_results[[paste(co,name)]]<<-z
 pz<-tryCatch(cox.zph(fit,transform='km')$table,error=function(e)NULL)
 if(!is.null(pz))ph[[paste(co,name)]]<<-data.frame(cohort=co,model=name,term=rownames(pz),chisq=pz[,1],df=pz[,2],p=pz[,3])
 return(fit)
}
for(co in names(cohorts)){
 d<-cohorts[[co]]
 cat('Analyzing',co,nrow(d),'\n')
 for(v in intersect(c('age','male','stage_num','grade_high','ALB','logAFP'),names(d))){
  x<-d[[v]];a<-x[!is.na(x)]
  descriptives[[paste(co,v)]]<-data.frame(cohort=co,variable=v,n_total=length(x),n_observed=length(a),n_missing=sum(is.na(x)),median=if(length(a))median(a)else NA,q1=if(length(a))quantile(a,.25)else NA,q3=if(length(a))quantile(a,.75)else NA,min=if(length(a))min(a)else NA,max=if(length(a))max(a)else NA)
 }
 for(v in intersect(c('OSARS_SD','High','age10','male','stage_advanced','grade_high','ALB','logAFP'),names(d)))fit_export(d,co,paste0('univariate_',v),v)
 adj<-if(grepl('CPTAC',co))c('OSARS_SD','age10','male')else c('OSARS_SD','age10','male','stage_advanced')
 fit_export(d,co,'multivariable_primary',adj)
 if(co=='TCGA-LIHC')fit_export(d,co,'multivariable_plus_grade',c(adj,'grade_high'))
 if(grepl('CPTAC',co))fit_export(d,co,'multivariable_plus_labs',c(adj,'ALB','logAFP'))
 for(v in intersect(c('High','High_train_cutoff'),names(d))){
  ss<-survdiff(as.formula(paste('Surv(time,event)~',v)),data=d)
  km_tests[[paste(co,v)]]<-data.frame(cohort=co,grouping=v,n=nrow(d),events=sum(d$event),logrank_chisq=ss$chisq,df=length(ss$n)-1,p=pchisq(ss$chisq,length(ss$n)-1,lower.tail=FALSE))
 }
 flin<-coxph(Surv(time,event)~OSARS_SD,data=d);fspl<-coxph(Surv(time,event)~splines::ns(OSARS_SD,df=3),data=d)
 lr<-2*(fspl$loglik[2]-flin$loglik[2]);nonlinear[[co]]<-data.frame(cohort=co,LR_chisq=lr,df=2,p=pchisq(lr,2,lower.tail=FALSE),test='three-df natural spline versus linear OSARS effect; exploratory diagnostic')
 write_out(d,paste0(co,'_clinical_analysis_rows.csv'))
}
write_out(do.call(rbind,all_results),'Cox_forest_long.csv')
write_out(do.call(rbind,ph),'Cox_PH_diagnostics.csv')
write_out(do.call(rbind,descriptives),'clinical_descriptive_statistics.csv')
write_out(do.call(rbind,nonlinear),'Cox_linearity_diagnostics.csv')
write_out(do.call(rbind,km_tests),'KM_logrank_tests.csv')
# Stage increment models are trained in TCGA ONLY and transferred to ICGC.
train<-cohorts[['TCGA-LIHC']]
train<-train[complete.cases(train[,c('stage_num','ROS_mean_rank','OSARS_frozen')]),]
forms<-list(Stage=Surv(time,event)~stage_num,Stage_plus_ROS=Surv(time,event)~stage_num+ROS_mean_rank,Stage_plus_OSARS=Surv(time,event)~stage_num+OSARS_frozen)
preds<-list();coefficients<-list()
for(n in names(forms)){
 m<-coxph(forms[[n]],data=train,ties='efron',x=TRUE,y=TRUE)
 coefficients[[n]]<-data.frame(model=n,term=names(coef(m)),coefficient=coef(m),n_training=nrow(train),training_events=sum(train$event))
 for(co in c('TCGA-LIHC','ICGC-LIRI')){
  dd<-cohorts[[co]];dd<-dd[complete.cases(dd[,c('stage_num','ROS_mean_rank','OSARS_frozen')]),]
  preds[[paste(n,co)]]<-data.frame(cohort=co,ID=dd$ID,time=dd$time,event=dd$event,predictor=n,score=predict(m,newdata=dd,type='lp',reference='zero'),role=if(co=='TCGA-LIHC')'apparent'else'TCGA-trained model, historically consulted external cohort')
 }
}
write_out(do.call(rbind,preds),'clinical_incremental_predictions.csv')
write_out(do.call(rbind,coefficients),'clinical_incremental_TCGA_coefficients.csv')
writeLines(capture.output(sessionInfo()),file.path(out,'logs/clinical_R_sessionInfo.txt'))
