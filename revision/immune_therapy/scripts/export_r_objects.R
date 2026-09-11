# Export existing objects only; do not refit models or alter original objects.
args <- commandArgs(trailingOnly=TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1]
module_dir <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), mustWork=FALSE)
root <- if (length(args)>=1) normalizePath(args[1],mustWork=FALSE) else module_dir
extra_library <- Sys.getenv("HCC_OS_R_LIBRARY")
if(nzchar(extra_library)).libPaths(c(extra_library,.libPaths()))
suppressPackageStartupMessages(library(gbm))
project <- if(length(args)>=3)normalizePath(args[3],mustWork=TRUE) else normalizePath(Sys.getenv("HCC_OS_LEGACY_PROJECT_ROOT",file.path(root,"../..")),mustWork=TRUE)
input_value <- if(length(args)>=2)args[2] else Sys.getenv("HCC_OS_IMMUNE_INPUT_DIR")
if(!nzchar(input_value))stop("Provide the immune source directory as argument 2 or HCC_OS_IMMUNE_INPUT_DIR.")
input <- normalizePath(input_value,mustWork=TRUE)
out <- file.path(root,"data/source")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(root,"qc"),recursive=TRUE,showWarnings=FALSE)
e <- new.env(); load(file.path(input,"ssgsea_result.Rdata"),envir=e)
write.csv(data.frame(Sample=colnames(e$re),t(e$re),check.names=FALSE),file.path(out,"ssgsea_patient_export.csv"),row.names=FALSE)
d <- new.env(); load(file.path(input,"GSE100797.Rdata"),envir=d)
write.csv(d$GSE100797_clin,file.path(out,"ACT_clinical_export.csv"),row.names=FALSE)
m <- new.env(); load(file.path(input,"StepCox_GBM_model_object.Rdata"),envir=m)
g <- m$gbm_model$selected_genes
write.csv(data.frame(Gene=g),file.path(out,"model_genes_export.csv"),row.names=FALSE)
norm <- toupper(trimws(rownames(d$GSE100797_expr)))
qc <- data.frame(N=ncol(d$GSE100797_expr),Features=length(g),Trees=m$gbm_model$best_iteration,
                 Duplicate_normalized_genes=sum(duplicated(norm)),NA_cells=sum(is.na(d$GSE100797_expr)),
                 Missing_genes=paste(setdiff(toupper(g),norm),collapse="|"))
write.csv(qc,file.path(out,"ACT_input_audit.csv"),row.names=FALSE)
print(qc)
t <- new.env();load(file.path(input,"QWEN0208_no_GSE14520.Rdata"),envir=t)
train <- t$list_train_vali_Data[["TCGA-LIHC"]][,g,drop=FALSE]
mu <- vapply(train,function(z) mean(as.numeric(z),na.rm=TRUE),numeric(1))
newdata <- as.data.frame(t(d$GSE100797_expr),check.names=FALSE)
colnames(newdata) <- toupper(trimws(colnames(newdata)))
newdata <- newdata[,!duplicated(colnames(newdata)),drop=FALSE]
for (z in setdiff(toupper(g),colnames(newdata))) newdata[[z]] <- NA_real_
newdata <- newdata[,toupper(g),drop=FALSE]
missing_fraction <- rowMeans(is.na(newdata))
for (z in g) newdata[[toupper(z)]][is.na(newdata[[toupper(z)]])] <- mu[[z]]
score <- as.numeric(predict(m$gbm_model$gbm_fit,newdata=newdata,n.trees=m$gbm_model$best_iteration,type="link"))
act <- data.frame(ID=rownames(newdata),Risk_Score=score,MissingFrac=missing_fraction,
                  Risk_Group=ifelse(score>median(score),"High","Low"))
act <- merge(act,d$GSE100797_clin,by="ID")
write.csv(act,file.path(out,"ACT_scores_current_model.csv"),row.names=FALSE)
tcga_score <- as.numeric(predict(m$gbm_model$gbm_fit,newdata=train,n.trees=m$gbm_model$best_iteration,type="link"))
tcga <- data.frame(Sample=rownames(train),Recomputed_risk=tcga_score)
arch <- read.csv(file.path(project,"FIG5_redone_drug_mutation/TCGA_LIHC_OSARS_groups.csv"))
check <- merge(tcga,arch,by="Sample")
stopifnot(nrow(check)==343,max(abs(check$Recomputed_risk-check$Risk_Score))<1e-10)
write.csv(data.frame(N=343,Max_abs_difference=max(abs(check$Recomputed_risk-check$Risk_Score)),GBM_version=as.character(packageVersion("gbm"))),file.path(out,"TCGA_reprediction_audit.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(root,"qc/R_session_info.txt"))
