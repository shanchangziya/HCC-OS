args<-commandArgs(trailingOnly=TRUE)
script_arg<-grep('^--file=',commandArgs(trailingOnly=FALSE),value=TRUE)[1]
out<-normalizePath(file.path(dirname(sub('^--file=','',script_arg)),'..'),mustWork=FALSE)
frozen_result<-if(length(args))args[1] else Sys.getenv('HCC_OS_FROZEN_RESULT')
if(!nzchar(frozen_result))stop('Provide res_no_GSE14520.Rdata as argument 1 or HCC_OS_FROZEN_RESULT.')
r<-new.env();load(frozen_result,envir=r)
res<-r$res;a<-res$ml.res[['GBM']];b<-res$ml.res[['StepCox[forward] + GBM']]
nm<-names(res$ml.res)
rows<-list()
add<-function(check,value)rows[[length(rows)+1]]<<-data.frame(check=check,value=as.character(value))
add('number_of_models',length(nm));add('single_models',sum(!grepl(' + ',nm,fixed=TRUE)))
add('two_stage_models',sum(grepl(' + ',nm,fixed=TRUE)))
add('GBM_and_StepCox_GBM_identical_feature_vector',identical(a$fit$var.names,b$fit$var.names))
add('GBM_and_StepCox_GBM_tree_counts',paste(a$best,b$best,sep=' / '))
for(co in names(res$riskscore[['GBM']])){
 x<-res$riskscore[['GBM']][[co]];y<-res$riskscore[['StepCox[forward] + GBM']][[co]]
 stopifnot(identical(x$ID,y$ID));delta<-max(abs(x$RS-y$RS))
 add(paste(co,'GBM_vs_StepCox_GBM_score_max_abs_difference'),delta);stopifnot(delta==0)
}
stopifnot(identical(a$fit$var.names,b$fit$var.names),a$best==b$best)
write.csv(do.call(rbind,rows),file.path(out,'data/processed/frozen_model_audit.csv'),row.names=FALSE)
