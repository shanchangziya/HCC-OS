suppressPackageStartupMessages({library(Seurat);library(UCell);library(AUCell);library(data.table)})
set.seed(20260905)
args <- commandArgs(trailingOnly = TRUE)
script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
src <- if (length(args) >= 1) args[[1]] else Sys.getenv("HCC_OS_DISCOVERY_ROOT")
validation_rds <- if (length(args) >= 2) args[[2]] else Sys.getenv("HCC_OS_GSE149614_RDS")
metadata_file <- if (length(args) >= 3) args[[3]] else Sys.getenv(
  "HCC_OS_GSE149614_METADATA",
  file.path(root, "data/raw/GSE149614_HCC.metadata.updated.txt.gz")
)
if (!nzchar(src) || !dir.exists(src)) stop("Set HCC_OS_DISCOVERY_ROOT or pass discovery root as argument 1")
if (!nzchar(validation_rds) || !file.exists(validation_rds)) stop("Set HCC_OS_GSE149614_RDS or pass it as argument 2")
if (!file.exists(metadata_file)) stop("GSE149614 metadata not found: ", metadata_file)
out <- file.path(root,"data/processed")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
ros <- readLines(file.path(root,"../model_validation/references/HALLMARK_ROS_msigdb_v7.0_genes.txt"))
ros <- unique(trimws(ros[nchar(ros)>0]))
mp <- read.csv(file.path(out,"discovery_NMF_genes.csv"))
programs <- split(mp$gene,mp$program)
programs$MP4_without_ROS <- setdiff(programs$MP4,ros)
stopifnot(length(ros)==49,length(programs$MP4)==58)
e <- new.env();load(file.path(src,"05osnmf.Rdata"),envir=e)
disc <- e$seu
f <- c("OS_AddModuleScore_1","OxStressOS_UCell_","OS_AUCell_OxStress")
dm <- disc@meta.data
dm$OS_composite_rederived <- rowMeans(scale(dm[,f]))
dm$OS_high_top20 <- dm$OS_score_mean_z>=quantile(dm$OS_score_mean_z,.8)
dm$OS_high_within_sample_top20 <- ave(dm$OS_score_mean_z,dm$orig.ident,FUN=function(x) x>=quantile(x,.8))>0
dm$cell <- rownames(dm)
dm$UMAP1 <- Embeddings(disc,"umap")[rownames(dm),1]
dm$UMAP2 <- Embeddings(disc,"umap")[rownames(dm),2]
dm$sample <- dm$orig.ident
dm$OS_composite <- dm$OS_score_mean_z
dm$OS_UCell <- dm$OxStressOS_UCell_
dm$OS_AddModuleScore <- dm$OS_AddModuleScore_1
dm$OS_AUCell <- dm$OS_AUCell_OxStress
write.csv(data.frame(metric=c("n_cells","n_tissue_samples","composite_rederived_max_abs_diff","top20_cutoff","n_top20","n_historical_valley"),value=c(nrow(dm),length(unique(dm$orig.ident)),max(abs(dm$OS_composite_rederived-dm$OS_score_mean_z)),unname(quantile(dm$OS_score_mean_z,.8)),sum(dm$OS_high_top20),sum(dm$OS_group_valley=="OS-high"))),file.path(out,"discovery_definition_audit.csv"),row.names=FALSE)
disc <- AddModuleScore_UCell(disc,features=programs,ncores=4,name="_recomputed",assay="RNA")
for(p in names(programs)) dm[[paste0(p,"_UCell")]] <- disc@meta.data[[paste0(p,"_recomputed")]]
write.csv(dm,file.path(out,"discovery_validation_ready.csv"),row.names=FALSE)
write.csv(data.frame(program=names(programs),genes=sapply(programs,length),overlap_ROS=sapply(programs,function(g)length(intersect(g,ros)))),file.path(out,"program_ROS_overlap.csv"),row.names=FALSE)
rm(e,disc);gc()

v <- readRDS(validation_rds)
meta <- fread(metadata_file,data.table=FALSE)
rownames(meta)<-paste0("GSE149614_",meta$Cell)
stopifnot(setequal(rownames(meta),colnames(v)))
meta <- meta[colnames(v),]
malignant_clusters<-c(3,4,12,15,17,19,22,24,27,42,43,47)
included <- meta$res.3 %in% malignant_clusters & meta$site=="Tumor"
v<-subset(v,cells=colnames(v)[included]);meta<-meta[included,]
v<-AddMetaData(v,meta)
v[["percent.mt"]]<-PercentageFeatureSet(v,pattern="^MT-")
v<-NormalizeData(v,normalization.method="LogNormalize",scale.factor=1e4,verbose=FALSE)
rg<-intersect(ros,rownames(v));stopifnot(length(rg)>=45)
v<-AddModuleScore(v,features=list(ROS=rg),name="ROS_AddModuleScore",assay="RNA",seed=20260905)
v<-AddModuleScore_UCell(v,features=c(list(ROS=rg),programs),ncores=4,name="_UCell",assay="RNA")
set.seed(20260905)
rankings<-AUCell_buildRankings(GetAssayData(v,assay="RNA",slot="data"),plotStats=FALSE,nCores=4,verbose=FALSE)
auc<-AUCell_calcAUC(list(ROS=rg),rankings,nCores=4,verbose=FALSE)
vm<-v@meta.data
vm$OS_AUCell<-as.numeric(getAUC(auc)["ROS",rownames(vm)])
vm$OS_AddModuleScore<-vm$ROS_AddModuleScore1
vm$OS_UCell<-vm$ROS_UCell
vm$OS_composite<-rowMeans(scale(vm[,c("OS_AddModuleScore","OS_UCell","OS_AUCell")]))
vm$OS_high_top20<-vm$OS_composite>=quantile(vm$OS_composite,.8)
vm$OS_high_within_sample_top20<-ave(vm$OS_composite,vm$patient,FUN=function(x)x>=quantile(x,.8))>0
vm$OS_UCell_high_within_sample_top20<-ave(vm$OS_UCell,vm$patient,FUN=function(x)x>=quantile(x,.8))>0
vm$technical_qc_pass<-vm$nFeature_RNA>=200 & vm$nFeature_RNA<=6000 & vm$percent.mt<20
vm$cell<-sub("^GSE149614_","",rownames(vm))
write.csv(vm,file.path(out,"validation_primary_scores.csv"),row.names=FALSE)
write.csv(data.frame(cohort="GSE149614",n_cells=nrow(vm),n_patients=length(unique(vm$patient)),n_ROS=length(rg),composite_cutoff=unname(quantile(vm$OS_composite,.8)),n_OS_high=sum(vm$OS_high_top20),n_QC_pass=sum(vm$technical_qc_pass)),file.path(out,"validation_primary_summary.csv"),row.names=FALSE)
saveRDS(v,file.path(out,"validation_primary_scored.rds"),compress=FALSE)
writeLines(capture.output(sessionInfo()),file.path(root,"audit/scoring_sessionInfo.txt"))
cat("FROZEN_PROGRAM_SCORING_COMPLETE\n")
