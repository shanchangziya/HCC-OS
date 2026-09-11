suppressPackageStartupMessages({library(Seurat);library(GeneNMF);library(UCell)})
set.seed(20260905)
script_arg<-commandArgs(trailingOnly=FALSE)
script_file<-sub("^--file=","",script_arg[grepl("^--file=",script_arg)][1])
root<-normalizePath(file.path(dirname(script_file),".."),mustWork=TRUE)
out<-file.path(root,"data/processed")
v<-readRDS(file.path(out,"validation_primary_scored.rds"))
DefaultAssay(v)<-"RNA"
v<-FindVariableFeatures(v,selection.method="vst",nfeatures=2000,verbose=FALSE)
objects<-SplitObject(v,split.by="patient")
# Retain the original discovery specification; do not select k on replication.
# GeneNMF::multiNMF internally calls set.seed(seed); its original default was 123.
# Make that actual historical specification explicit, without changing the fit.
p<-multiNMF(objects,assay="RNA",k=4:9,min.exp=.05,seed=123)
m<-getMetaPrograms(p,metric="cosine",specificity.weight=3,weight.explained=.3,nMP=4)
saveRDS(list(programs=p,metaprograms=m),file.path(out,"validation_de_novo_NMF.rds"))
genes<-m$metaprograms.genes
write.csv(do.call(rbind,lapply(names(genes),function(n)data.frame(program=n,gene=genes[[n]]))),file.path(out,"validation_de_novo_NMF_genes.csv"),row.names=FALSE)
write.csv(m$metaprograms.metrics,file.path(out,"validation_de_novo_NMF_metrics.csv"))
d<-read.csv(file.path(out,"discovery_NMF_genes.csv"));dg<-split(d$gene,d$program)
# Condition the overlap null on the genes actually entered into the saved NMF W.
# Preserve measured-gene Jaccard separately as a descriptive similarity statistic.
source(file.path(root,"scripts/06_recompute_NMF_overlap.R"))
recompute_NMF_overlap(p,genes,dg,rownames(v),out,file.path(root,"audit"),123L)
v<-AddModuleScore_UCell(v,features=genes,ncores=4,name="_deNovo",assay="RNA")
scores<-v@meta.data[,c("Cell","patient",paste0(names(genes),"_deNovo"))]
write.csv(scores,file.path(out,"validation_de_novo_NMF_scores.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(root,"audit/NMF_sessionInfo.txt"))
cat("VALIDATION_NMF_COMPLETE\n")
