suppressPackageStartupMessages(library(Seurat))
args <- commandArgs(trailingOnly = TRUE)
script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
module_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
source_dir <- if (length(args) >= 1) args[[1]] else Sys.getenv("HCC_OS_DISCOVERY_ROOT")
out <- if (length(args) >= 2) args[[2]] else module_dir
validation_rds <- if (length(args) >= 3) args[[3]] else Sys.getenv("HCC_OS_GSE149614_RDS")
if (!nzchar(source_dir) || !dir.exists(source_dir)) {
  stop("Provide the discovery pipeline directory as argument 1 or HCC_OS_DISCOVERY_ROOT")
}
if (!nzchar(validation_rds) || !file.exists(validation_rds)) {
  stop("Provide the GSE149614 Seurat RDS as argument 3 or HCC_OS_GSE149614_RDS")
}
dir.create(file.path(out, "data/processed"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out, "audit"), recursive=TRUE, showWarnings=FALSE)
probe <- function(obj, label) {
  cat("\nOBJECT", label, "class", class(obj), "dimension", dim(obj), "\n")
  cat("ASSAYS", names(obj@assays), "REDUCTIONS", names(obj@reductions), "\n")
  print(colnames(obj@meta.data))
  print(utils::head(obj@meta.data, 2))
  for (c in intersect(c("orig.ident", "sample", "patient", "Patient", "dataset", "group", "celltype", "cell_type", "malignant_status", "OS_group", "OS_group_valley"), colnames(obj@meta.data))) {
    cat("METADATA", c, "\n"); print(table(obj@meta.data[[c]],useNA="ifany"))
  }
  write.csv(cbind(cell=colnames(obj),obj@meta.data), file.path(out,"data/processed",paste0(label,"_metadata_original.csv")), row.names=FALSE)
  cat("GENE HEAD", head(rownames(obj)), "\nCELL HEAD", head(colnames(obj)), "\n")
}
v <- readRDS(validation_rds)
probe(v,"GSE149614")
rm(v);gc()
e <- new.env(); print(load(file.path(source_dir,"05osnmf.Rdata"),envir=e))
for(n in ls(e)) if(inherits(e[[n]],"Seurat")) probe(e[[n]],paste0("discovery_",n))
if(exists("geneNMF.metaprograms",e)) {
 m <- e$geneNMF.metaprograms
 cat("NMF FIELDS", names(m), "\n"); print(m$metaprograms.metrics)
 genes <- m$metaprograms.genes; print(lapply(genes,head))
 write.csv(do.call(rbind,lapply(names(genes), function(n) data.frame(program=n,gene=genes[[n]]))),file.path(out,"data/processed/discovery_NMF_genes.csv"),row.names=FALSE)
}
print(sessionInfo())
