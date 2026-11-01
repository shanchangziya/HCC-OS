# Export only the numeric fields needed for inference. Use the original trusted
# RData file containing gbm_model; no fitting or patient data export is required.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript export_frozen_model.R ORIGINAL.Rdata NEW_MODEL_DIRECTORY")
}

export_model <- function(source, destination) {
  if (file.exists(destination)) stop("Choose a new output directory.")
  env <- new.env(parent = baseenv())
  load(source, envir = env)
  model <- env$gbm_model
  fit <- model$gbm_fit
  genes <- fit$var.names
  n <- length(fit$trees)
  stopifnot(
    inherits(fit, "gbm"), identical(as.character(model$selected_genes), genes),
    length(genes) == 128L, !anyDuplicated(genes), n == 1696L,
    fit$n.trees == n, model$best_iteration == n, fit$initF == 0,
    identical(fit$distribution$name, "coxph"), all(fit$var.type == 0),
    length(fit$c.splits) == 0L
  )

  stage <- tempfile("osars-export-", tmpdir = dirname(destination))
  if (!dir.create(stage)) stop("Could not create the temporary output directory.")
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  writeLines(c("format\tosars-numeric-gbm-v1", "n_features\t128",
               "n_trees\t1696", "initF\t0", "prediction_scale\tlink"),
             file.path(stage, "model.tsv"))
  writeLines(genes, file.path(stage, "genes.txt"))

  path <- file.path(stage, "trees.tsv")
  # R gbm uses zero-based child and split-variable indices; keep them verbatim.
  writeLines("tree\tnode\tfeature\tsplit_or_leaf\tleft\tright\tmissing", path)
  for (i in seq_len(n)) {
    tr <- fit$trees[[i]]
    terminal <- tr[[1L]] == -1L
    stopifnot(all(tr[[2L]][terminal] == tr[[8L]][terminal]),
              all(is.finite(tr[[2L]])))
    lines <- sprintf("%d\t%d\t%d\t%.17g\t%d\t%d\t%d", i,
                     seq_along(tr[[1L]]) - 1L, tr[[1L]], tr[[2L]],
                     tr[[3L]], tr[[4L]], tr[[5L]])
    cat(lines, file = path, sep = "\n", append = TRUE)
  }
  if (!file.rename(stage, destination)) stop("Could not publish the output directory.")
}
export_model(args[[1L]], args[[2L]])
