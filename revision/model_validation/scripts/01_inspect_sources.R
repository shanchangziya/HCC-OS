args<-commandArgs(trailingOnly=TRUE)
script_arg<-grep('^--file=',commandArgs(trailingOnly=FALSE),value=TRUE)[1]
module_dir<-normalizePath(file.path(dirname(sub('^--file=','',script_arg)),'..'),mustWork=FALSE)
root <- if(length(args)>=1) args[1] else Sys.getenv('HCC_OS_DISCOVERY_ROOT')
out <- if(length(args)>=2) args[2] else module_dir
if(!nzchar(root))stop('Provide the archived discovery-analysis root as argument 1 or HCC_OS_DISCOVERY_ROOT.')
dir.create(file.path(out,'data/processed'), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out,'data/raw'), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out,'logs'), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out,'references'), recursive=TRUE, showWarnings=FALSE)
sink(file.path(out,'logs/source_structure.txt'))
for(f in c('QWEN0208.Rdata','QWEN0208_no_GSE14520.Rdata','StepCox_GBM_model_object.Rdata','bulk/exp1surv1.rdata','bulk/exp2surv2.rdata','bulk/TCGA-LIHC.Rdata','临床治疗/exp1_surv1.Rdata')) {
 if(!file.exists(file.path(root,f))) next
 e <- new.env(); load(file.path(root,f),envir=e)
 cat('\nFILE: ',f,'\n')
 for(n in ls(e)) {x<-get(n,e);cat(n, class(x), dim(x),'\n'); if(is.list(x))str(x,max.level=1,list.len=6);if(is.matrix(x)) {print(x[seq_len(min(3,nrow(x))),seq_len(min(6,ncol(x)))]);print(quantile(x,na.rm=T))}}
}
for(f in c('TCGA临床信息.txt','bulk/TCGA临床信息.txt')) {if(!file.exists(file.path(root,f)))next;cat('\n',f,'\n'); x<-read.delim(file.path(root,f)); print(names(x)); print(head(x,2))}
e<-new.env();load(file.path(root,'res_no_GSE14520.Rdata'),envir=e)
cat('\nRES\n');str(e$res,max.level=2,list.len=4)
m<-e$res$ml.res[['StepCox[forward] + GBM']];cat('\nMODEL\n');str(m,max.level=2,list.len=15)
cat('\nRISK\n');str(e$res$riskscore[['StepCox[forward] + GBM']],max.level=2)
write.csv(e$res$Cindex.res,file.path(out,'data/processed/original_117_model_performance.csv'),row.names=FALSE)
saveRDS(e$res$riskscore[['StepCox[forward] + GBM']],file.path(out,'data/raw/frozen_risks.rds'))
if(requireNamespace('Mime1',quietly=TRUE)){writeLines(capture.output(Mime1::ML.Dev.Prog.Sig),file.path(out,'references/Mime1_ML.Dev.Prog.Sig_source.txt'))}
print(sessionInfo())
sink()
