#!/usr/bin/env Rscript
# Numerical backend only: no graphics calls or plot devices.
args <- commandArgs(trailingOnly=TRUE)
out <- normalizePath(args[1], mustWork=TRUE)
B <- as.integer(args[2])
if(length(args)>=3) .libPaths(c(args[3], .libPaths()))
suppressPackageStartupMessages({library(survival); library(glmnet); library(jsonlite)})
write_table <- function(x, name) write.csv(x,file.path(out,paste0(name,".csv")),row.names=FALSE,na="")
write_json_safe <- function(x,name) write_json(x,file.path(out,paste0(name,".json")),
  pretty=TRUE, auto_unbox=TRUE, na="null", digits=15)
dat <- read.csv(file.path(out,"analysis_input.csv"),check.names=FALSE)
features <- paste0("resnet",0:2047)
stopifnot(all(features %in% names(dat)), !anyDuplicated(dat$patient_id),
          sum(dat$split=="Train")==230, sum(dat$split=="Test")==100,
          all(dat$time_days>0),all(dat$event %in% c(0,1)))
train_idx <- which(dat$split=="Train")
set.seed(2026)
legacy_seed_match <- setequal(sample(seq_len(nrow(dat)),floor(.7*nrow(dat))),train_idx)
train <- dat[train_idx,,drop=FALSE]
test_idx <- which(dat$split=="Test")
X <- as.matrix(dat[,features])
storage.mode(X)<-"double"
screen_features <- function(rows) {
  target <- dat$molecular_OSARS[rows]
  tt <- lapply(features,function(f) {
    x <- X[rows,f]
    if(sd(x)==0 || sd(target)==0) return(data.frame(feature=f,n=length(x),r=NA,p=NA))
    ct<-suppressWarnings(cor.test(x,target,method="pearson"))
    data.frame(feature=f,n=length(x),r=unname(ct$estimate),p=ct$p.value)
  })
  tab <- do.call(rbind,tt)
  tab$fdr <- p.adjust(tab$p,method="BH")
  tab$selected <- is.finite(tab$r) & is.finite(tab$p) & abs(tab$r)>.2 & tab$p<.05
  tab
}
# Fold-specific screening and standardization prevent validation-fold information
# from entering training-fold transformations during lambda selection.
set.seed(20260905)
foldid <- integer(nrow(train))
for(e in c(0,1)) {
  ii<-which(train$event==e)
  foldid[ii]<-sample(rep(1:10,length.out=length(ii)))
}
write_table(data.frame(patient_id=train$patient_id,event=train$event,fold=foldid),"training_cv_folds")
ratios <- exp(seq(log(1),log(.001),length.out=60))
fold_loss <- matrix(NA_real_,10,length(ratios))
fold_meta<-list(); fold_features<-list()
# Held-out Breslow partial likelihood, divided by validation-fold event count.
# Additive constants irrelevant to lambda choice are omitted.
partial_deviance <- function(lp,time,event) {
  if(sum(event)==0) return(NA_real_)
  lp<-lp-max(lp)
  ll<-sum(vapply(which(event==1),function(i)
    lp[i]-log(sum(exp(lp[time>=time[i]]))),numeric(1)))
  -2*ll/sum(event)
}
for(k in 1:10) {
  local_train<-which(foldid!=k); local_valid<-which(foldid==k)
  sel<-screen_features(train_idx[local_train])
  ff<-sel$feature[sel$selected]
  fold_features[[k]]<-data.frame(fold=k,sel)
  yy<-Surv(train$time_days[local_train],train$event[local_train])
  if(length(ff)>=2) {
    probe<-glmnet(X[train_idx[local_train],ff,drop=FALSE],yy,family="cox",alpha=1,
                   standardize=TRUE,nlambda=30,lambda.min.ratio=.001,maxit=100000)
    lambda_max<-max(probe$lambda)
    model<-glmnet(X[train_idx[local_train],ff,drop=FALSE],yy,family="cox",alpha=1,
                   standardize=TRUE,lambda=lambda_max*ratios,maxit=100000)
    lp<-as.matrix(predict(model,newx=X[train_idx[local_valid],ff,drop=FALSE],
                           s=lambda_max*ratios,type="link"))
  } else {
    lambda_max<-NA_real_
    lp<-matrix(0,length(local_valid),length(ratios))
  }
  fold_loss[k,]<-apply(lp,2,partial_deviance,time=train$time_days[local_valid],
                     event=train$event[local_valid])
  fold_meta[[k]]<-data.frame(fold=k,n_train=length(local_train),n_validation=length(local_valid),
    events_validation=sum(train$event[local_valid]),selected_features=length(ff),lambda_max=lambda_max)
  cat("CV fold",k,"features",length(ff),"completed\n")
}
write_table(do.call(rbind,fold_meta),"cv_fold_summary")
write_table(do.call(rbind,fold_features),"cv_fold_feature_screening")
cvm<-colMeans(fold_loss)
cvse<-apply(fold_loss,2,sd)/sqrt(nrow(fold_loss))
best<-which.min(cvm)
one_se<-min(which(cvm<=cvm[best]+cvse[best]))
screen<-screen_features(train_idx)
write_table(screen,"train_feature_correlations")
selected<-screen$feature[screen$selected]
if(length(selected)<2) stop("Insufficient training-only correlated features: report without model.")
train_y<-Surv(train$time_days,train$event)
probe<-glmnet(X[train_idx,selected,drop=FALSE],train_y,family="cox",alpha=1,
              standardize=TRUE,nlambda=30,lambda.min.ratio=.001,maxit=100000)
lambda_max<-max(probe$lambda)
lambda<-lambda_max*ratios[best]
path<-glmnet(X[train_idx,selected,drop=FALSE],train_y,family="cox",alpha=1,
             standardize=TRUE,lambda=lambda_max*ratios,maxit=100000)
coeff<-as.matrix(coef(path,s=lambda))[,1]
train_lp<-as.numeric(predict(path,newx=X[train_idx,selected,drop=FALSE],s=lambda,type="link"))
stopifnot(isTRUE(all.equal(as.numeric(X[train_idx,selected,drop=FALSE]%*%coeff),train_lp,tolerance=1e-8)))
cutoff<-median(train_lp)
if(sd(train_lp)==0) stop("Prespecified CV selected null model: no prognostic discrimination to plot.")
coef_table<-data.frame(feature=selected,coefficient_raw_scale=coeff,
  train_mean=colMeans(X[train_idx,selected,drop=FALSE]),
  train_sd=apply(X[train_idx,selected,drop=FALSE],2,sd),nonzero=coeff!=0)
coef_table$coefficient_per_training_sd<-coef_table$coefficient_raw_scale*coef_table$train_sd
write_table(coef_table,"lasso_coefficients")
write_table(data.frame(lambda=lambda_max*ratios,lambda_fraction=ratios,
  log_lambda=log(lambda_max*ratios),cv_mean_deviance=cvm,cv_se=cvse,
  cv_lower=cvm-cvse,cv_upper=cvm+cvse,selected_min=seq_along(ratios)==best,
  selected_1se=seq_along(ratios)==one_se),"cv_curve")
coef_path<-as.matrix(coef(path,s=lambda_max*ratios))
write_table(data.frame(feature=rep(rownames(coef_path),ncol(coef_path)),
  lambda=rep(lambda_max*ratios,each=nrow(coef_path)),
  log_lambda=rep(log(lambda_max*ratios),each=nrow(coef_path)),
  coefficient=as.vector(coef_path)),"coefficient_path")
train$pathology_score<-train_lp
clinical_vars<-c("stage_advanced","grade_high","sex_male")
clinical_models<-list()
clinical_coefs<-list()
if(all(clinical_vars %in% names(train))) {
  complete<-complete.cases(train[,clinical_vars])
  if(sum(train$event[complete])>=30) {
    clinical_models$Clinical<-coxph(Surv(time_days,event)~stage_advanced+grade_high+sex_male,
      data=train[complete,],x=TRUE,y=TRUE,model=TRUE)
    clinical_models$Clinical_plus_pathology<-coxph(Surv(time_days,event)~pathology_score+
      stage_advanced+grade_high+sex_male,data=train[complete,],x=TRUE,y=TRUE,model=TRUE)
    for(nm in names(clinical_models)) {
      s<-summary(clinical_models[[nm]])
      clinical_coefs[[nm]]<-data.frame(model=nm,term=rownames(s$coefficients),
        coefficient=s$coefficients[,"coef"],hr=s$conf.int[,"exp(coef)"],
        lower=s$conf.int[,"lower .95"],upper=s$conf.int[,"upper .95"],
        p=s$coefficients[,"Pr(>|z|)"],training_n=sum(complete),training_events=sum(train$event[complete]))
    }
  }
}
if(length(clinical_coefs)) write_table(do.call(rbind,clinical_coefs),"clinical_model_coefficients")
locked<-list(frozen_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
  split_source="pathdat.Rdata dat_all$Set",legacy_seed_2026_reproduces_membership=legacy_seed_match,
  train_n=nrow(train),screening_n_features=length(selected),nonzero_n=sum(coeff!=0),
  screening="Training-only Pearson |r|>0.2 and nominal P<0.05; BH FDR reported",
  cv="10 stratified folds; fold-specific correlation screening and glmnet standardization; minimum mean held-out Breslow deviance per event",
  cv_seed=20260905,lambda_fraction=ratios[best],lambda=lambda,lambda_1se=lambda_max*ratios[one_se],
  lambda_max=lambda_max,training_median_cutoff=cutoff,group_rule="High if score > training median; else Low",
  coefficient_scale="Original raw ResNet feature scale; glmnet standardization learned on training only",
  clinical_models=names(clinical_models),clinical_variables=clinical_vars,
  auc_horizon_days=c(365.25,1095.75,1826.25),bootstrap_replicates=B,
  bootstrap_seed=20260906,limitations=c("Historical test cohort previously examined; corrected internal reassessment",
    "Molecular OSARS target was itself developed on TCGA; target non-independence remains",
    "No external pathology cohort; nominal exploratory correlation screening",
    "Bootstrap CI is conditional on fixed fitted model and does not capture training/model-selection uncertainty"))
write_json_safe(locked,"frozen_model_specification")
saveRDS(list(model=path,lambda=lambda,features=selected,coefficients=coeff,cutoff=cutoff,
  clinical_models=clinical_models,train_ids=train$patient_id,specification=locked),
  file.path(out,"frozen_pathology_model.rds"))
cat("MODEL FROZEN; beginning fixed test evaluation\n")
# No model, hyperparameter, feature, or cutoff changes are made after this line.
dat$pathology_score<-as.numeric(predict(path,newx=X[,selected,drop=FALSE],s=lambda,type="link"))
dat$risk_group<-ifelse(dat$pathology_score>cutoff,"High","Low")
dat$risk_group<-factor(dat$risk_group,levels=c("Low","High"))
for(nm in names(clinical_models)) dat[[paste0("score_",nm)]]<-
  as.numeric(predict(clinical_models[[nm]],newdata=dat,type="lp",reference="zero"))
export_cols<-c("patient_id","ID","split","time_days","event","n_tiles","molecular_OSARS", "pathology_score",
  "risk_group",intersect(c("Stage_Group","Grade_Group","Sex",clinical_vars),names(dat)),
  paste0("score_",names(clinical_models)))
write_table(dat[,export_cols],"frozen_patient_predictions")
write_table(aggregate(event~split+risk_group,dat,function(x)c(n=length(x),events=sum(x))),"cohort_group_counts")
km_all<-list(); censor_all<-list(); risk_all<-list(); group_results<-list(); ph_all<-list()
risk_times<-seq(0,5*365.25,by=365.25)
for(sp in c("Train","Test")) {
  d<-dat[dat$split==sp,]
  sf<-survfit(Surv(time_days,event)~risk_group,data=d,conf.type="log-log")
  su<-summary(sf,censored=TRUE)
  km<-data.frame(split=sp,group=sub("risk_group=","",as.character(su$strata)),
    time_days=su$time,survival=su$surv,lower=su$lower,upper=su$upper,
    n_risk=su$n.risk,n_event=su$n.event,n_censor=su$n.censor)
  start<-data.frame(split=sp,group=c("Low","High"),time_days=0,survival=1,lower=1,upper=1,
    n_risk=as.integer(table(d$risk_group)),n_event=0,n_censor=0)
  km<-rbind(start,km); km<-km[order(km$group,km$time_days),]
  km$time_years<-km$time_days/365.25
  km_all[[sp]]<-km
  censor_all[[sp]]<-km[km$n_censor>0,c("split","group","time_days","time_years","survival","n_censor")]
  rs<-summary(sf,times=risk_times,extend=TRUE)
  risk_all[[sp]]<-data.frame(split=sp,group=sub("risk_group=","",as.character(rs$strata)),
    time_days=rs$time,time_years=rs$time/365.25,n_risk=rs$n.risk,n_event_interval=rs$n.event,n_censor_interval=rs$n.censor)
  lr<-survdiff(Surv(time_days,event)~risk_group,data=d)
  fit<-coxph(Surv(time_days,event)~risk_group,data=d,x=TRUE,y=TRUE)
  ss<-summary(fit)
  group_results[[sp]]<-data.frame(split=sp,n=nrow(d),events=sum(d$event),n_low=sum(d$risk_group=="Low"),
    n_high=sum(d$risk_group=="High"),cutoff=cutoff,logrank_chisq=lr$chisq,logrank_p=pchisq(lr$chisq,1,lower.tail=FALSE),
    hr_high_vs_low=ss$conf.int[1,"exp(coef)"],hr_lower=ss$conf.int[1,"lower .95"],
    hr_upper=ss$conf.int[1,"upper .95"],cox_p=ss$coefficients[1,"Pr(>|z|)"])
  for(term in c("risk_group","pathology_score")) {
    f<-coxph(as.formula(paste("Surv(time_days,event)~",term)),data=d,x=TRUE,y=TRUE)
    zz<-tryCatch(cox.zph(f)$table,error=function(e)NULL)
    if(!is.null(zz)) ph_all[[paste(sp,term)]]<-data.frame(split=sp,model=term,term=rownames(zz),zz)
  }
}
write_table(do.call(rbind,km_all),"km_coordinates")
write_table(do.call(rbind,censor_all),"km_censor_coordinates")
write_table(do.call(rbind,risk_all),"km_risk_table")
write_table(do.call(rbind,group_results),"survival_group_statistics")
write_table(do.call(rbind,ph_all),"proportional_hazards_diagnostics")

times<-c(365.25,1095.75,1826.25)
has_timeroc<-requireNamespace("timeROC",quietly=TRUE)
if(!has_timeroc) cat("timeROC unavailable; equivalent marginal-KM IPCW calculation used.\n")
# Cumulative/dynamic AUC: cases T<=t, event=1; controls T>t. Censored
# observations before t are excluded. Case weights use G(T_i-); controls G(t).
km_ipcw <- function(time,event,score,horizon,return_curve=FALSE) {
  case<-which(time<=horizon & event==1); control<-which(time>horizon)
  if(length(case)<2 || length(control)<2) return(list(auc=NA_real_))
  censorfit<-survfit(Surv(time,1-event)~1)
  Gbefore<-function(x) {
    ii<-which(censorfit$time<x)
    if(length(ii)) tail(censorfit$surv[ii],1) else 1
  }
  ww<-1/vapply(time[case],Gbefore,numeric(1))
  if(any(!is.finite(ww))) return(list(auc=NA_real_))
  delta<-outer(score[case],score[control],"-")
  auc<-sum((delta>0)*ww + .5*(delta==0)*ww)/(sum(ww)*length(control))
  curve<-NULL
  if(return_curve) {
    thresholds<-c(Inf,sort(unique(score),decreasing=TRUE),-Inf)
    curve<-data.frame(threshold=thresholds,
      sensitivity=vapply(thresholds,function(cut)sum(ww[score[case]>=cut])/sum(ww),numeric(1)),
      false_positive_rate=vapply(thresholds,function(cut)mean(score[control]>=cut),numeric(1)))
  }
  list(auc=auc,curve=curve)
}
get_cindex<-function(time,event,score) {
  tryCatch(unname(concordance(Surv(time,event)~score,reverse=TRUE)$concordance),error=function(e)NA_real_)
}
get_auc<-function(time,event,score) {
  vapply(times,function(t) km_ipcw(time,event,score,t)$auc,numeric(1))
}
metrics<-list(); roc_data<-list(); bootstrap_tables<-list(); paired_delta<-list(); method_checks<-list()
set.seed(20260906)
for(sp in c("Train","Test")) {
  base<-dat[dat$split==sp,]
  subsets<-list(All=base)
  if(length(clinical_models)) subsets$Clinical_complete<-base[complete.cases(base[,clinical_vars]),]
  for(subset_name in names(subsets)) {
    d<-subsets[[subset_name]]
    models<-c(Pathology="pathology_score")
    if(subset_name=="Clinical_complete") models<-c(models,setNames(paste0("score_",names(clinical_models)),names(clinical_models)))
    bootstrap_idx<-replicate(B,sample.int(nrow(d),nrow(d),replace=TRUE),simplify=FALSE)
    boot_model<-list()
    for(model_name in names(models)) {
      score<-d[[models[[model_name]]]]
      cc<-concordance(Surv(time_days,event)~score,data=d,reverse=TRUE)
      point<-c(C_index=unname(cc$concordance),setNames(get_auc(d$time_days,d$event,score),c("AUC_1y","AUC_3y","AUC_5y")))
      boots<-t(vapply(bootstrap_idx,function(ii) c(get_cindex(d$time_days[ii],d$event[ii],score[ii]),
        get_auc(d$time_days[ii],d$event[ii],score[ii])),numeric(4)))
      colnames(boots)<-names(point)
      boot_model[[model_name]]<-boots
      boot_df<-data.frame(split=sp,subset=subset_name,model=model_name,replicate=seq_len(B),boots)
      bootstrap_tables[[paste(sp,subset_name,model_name)]]<-boot_df
      for(m in names(point)) {
        valid<-boots[is.finite(boots[,m]),m]
        ci<-if(length(valid)>=max(100,.8*B)) quantile(valid,c(.025,.975),names=FALSE) else c(NA,NA)
        horizon<-if(m=="C_index") NA_real_ else times[match(m,c("AUC_1y","AUC_3y","AUC_5y"))]
        metrics[[paste(sp,subset_name,model_name,m)]]<-data.frame(split=sp,subset=subset_name,model=model_name,
          metric=m,estimate=point[[m]],lower=ci[1],upper=ci[2],n=nrow(d),events=sum(d$event),
          bootstrap_valid=length(valid),bootstrap_requested=B,horizon_days=horizon,
          horizon_cases=if(is.na(horizon))NA else sum(d$time_days<=horizon & d$event==1),
          horizon_controls=if(is.na(horizon))NA else sum(d$time_days>horizon),
          cindex_influence_se=if(m=="C_index")sqrt(cc$var) else NA,
          ci_method="paired patient bootstrap percentile; frozen predictions; marginal KM IPCW for AUC")
      }
      for(j in seq_along(times)) {
        rr<-km_ipcw(d$time_days,d$event,score,times[j],TRUE)
        if(!is.null(rr$curve)) roc_data[[paste(sp,subset_name,model_name,j)]]<-
          data.frame(split=sp,subset=subset_name,model=model_name,time_years=c(1,3,5)[j],rr$curve)
      }
      if(has_timeroc) {
        tr<-tryCatch(timeROC::timeROC(T=d$time_days,delta=d$event,marker=score,cause=1,
          weighting="marginal",times=times,iid=FALSE),error=function(e)NULL)
        if(!is.null(tr)) method_checks[[paste(sp,subset_name,model_name)]]<-data.frame(split=sp,
          subset=subset_name,model=model_name,time_years=c(1,3,5),our_auc=point[2:4],timeROC_auc=tr$AUC,
          abs_difference=abs(point[2:4]-tr$AUC))
      }
      cat(sp,subset_name,model_name,"C-index",point[1],"AUC",point[2:4],"\n")
    }
    if(all(c("Clinical","Clinical_plus_pathology") %in% names(boot_model))) {
      delta<-boot_model$Clinical_plus_pathology-boot_model$Clinical
      for(m in colnames(delta)) {
        vals<-delta[is.finite(delta[,m]),m]
        pc<-metrics[[paste(sp,subset_name,"Clinical_plus_pathology",m)]]$estimate-
            metrics[[paste(sp,subset_name,"Clinical",m)]]$estimate
        ci<-quantile(vals,c(.025,.975),names=FALSE)
        paired_delta[[paste(sp,subset_name,m)]]<-data.frame(split=sp,subset=subset_name,
          comparison="Clinical plus pathology minus clinical",metric=m,difference=pc,
          lower=ci[1],upper=ci[2],bootstrap_valid=length(vals))
      }
    }
  }
}
write_table(do.call(rbind,metrics),"performance_metrics")
write_table(do.call(rbind,roc_data),"roc_coordinates")
write_table(do.call(rbind,bootstrap_tables),"bootstrap_performance")
if(length(paired_delta))write_table(do.call(rbind,paired_delta),"paired_incremental_performance")
if(length(method_checks))write_table(do.call(rbind,method_checks),"auc_implementation_crosscheck")
write_table(data.frame(package=c("R","survival","glmnet","jsonlite","timeROC"),
  version=c(R.version.string,as.character(packageVersion("survival")),as.character(packageVersion("glmnet")),
  as.character(packageVersion("jsonlite")),if(has_timeroc)as.character(packageVersion("timeROC"))else NA)),"software_versions")
capture.output(sessionInfo(),file=file.path(out,"R_session_info.txt"))
cat("Completed numerical exports without generating plots.\n")
