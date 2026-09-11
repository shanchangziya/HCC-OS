args <- commandArgs(trailingOnly=TRUE)
src <- normalizePath(args[1], mustWork=TRUE)
out <- args[2]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
sink(file.path(out, "source_object_inventory.txt"))
cat("UTC:", format(Sys.time(), tz="UTC"), "\n")
for (f in c("Resnet.Rdata", "pathdat.Rdata")) {
  cat("\nSOURCE", f, "\n")
  env <- new.env(parent=globalenv())
  load(file.path(src, f), envir=env)
  for (nm in ls(env)) {
    x <- env[[nm]]
    cat("OBJECT", nm, "CLASS", paste(class(x),collapse=","),
        "DIM", paste(dim(x),collapse="x"), "LENGTH", length(x), "\n")
    if (is.data.frame(x) || is.matrix(x)) {
      cat("COLNAMES", paste(head(colnames(x), 160), collapse=" | "), "\n")
      cat("ROWNAMES", paste(head(rownames(x), 10), collapse=" | "), "\n")
      print(x[seq_len(min(3,nrow(x))),seq_len(min(10,ncol(x))),drop=FALSE])
    } else if (is.atomic(x) && length(x) <= 350) {
      print(x)
    } else if (is.list(x)) {
      cat("NAMES", paste(names(x),collapse=" | "), "\n")
      str(x, max.level=1, vec.len=4, list.len=50)
    }
  }
}
for (p in c("survival", "glmnet", "timeROC", "jsonlite", "caret")) {
  cat("PACKAGE", p, requireNamespace(p, quietly=TRUE), "\n")
  if(requireNamespace(p, quietly=TRUE)) cat(as.character(packageVersion(p)), "\n")
}
sessionInfo()
sink()
