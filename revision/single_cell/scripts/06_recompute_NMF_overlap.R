# Correct the overlap null universe using an already fitted NMF model.
# This script never fits or refits NMF and never draws a plot.

recompute_NMF_overlap <- function(programs, validation_genes, discovery_genes,
                                 measured_universe, out, audit,
                                 actual_NMF_seed=123L) {
  dir.create(audit, recursive=TRUE, showWarnings=FALSE)
  eligible_sets <- lapply(programs, function(model) rownames(model$w))
  stopifnot(length(eligible_sets)>0, all(vapply(eligible_sets,length,integer(1))>0))
  universe <- eligible_sets[[1]]
  stopifnot(!anyDuplicated(universe), all(vapply(eligible_sets,
    function(x)setequal(x,universe),logical(1))))
  stopifnot(all(unlist(validation_genes) %in% universe))

  # Keep the originally exported full-measured-gene Fisher results as audit only.
  prior <- file.path(out,"cross_cohort_NMF_overlap.csv")
  legacy <- file.path(audit,"cross_cohort_NMF_overlap_full_measured_gene_background_original.csv")
  if(file.exists(prior) && !file.exists(legacy)) {
    old <- read.csv(prior)
    if(!"background_n" %in% names(old)) file.copy(prior,legacy,overwrite=FALSE)
  }

  tab <- do.call(rbind,lapply(names(discovery_genes),function(a) {
    do.call(rbind,lapply(names(validation_genes),function(b) {
      # Raw Jaccard remains the original measured-gene descriptive statistic.
      raw_x <- intersect(discovery_genes[[a]],measured_universe)
      raw_y <- intersect(validation_genes[[b]],measured_universe)
      raw_k <- length(intersect(raw_x,raw_y))
      x <- intersect(discovery_genes[[a]],universe)
      y <- intersect(validation_genes[[b]],universe)
      k <- length(intersect(x,y))
      contingency <- matrix(c(k,length(x)-k,length(y)-k,
        length(universe)-length(union(x,y))),nrow=2)
      data.frame(discovery=a,validation=b,
        n_discovery=length(x),n_validation=length(y),overlap=k,
        jaccard=raw_k/length(union(raw_x,raw_y)),
        P=fisher.test(contingency,alternative="greater")$p.value,
        background_n=length(universe),
        background_definition="Genes actually eligible for validation NMF, from all saved W row names",
        n_discovery_raw_measured=length(raw_x),n_validation_raw_measured=length(raw_y),
        overlap_raw_measured=raw_k,jaccard_eligible=k/length(union(x,y)),
        jaccard_unrestricted_gene_sets=length(intersect(discovery_genes[[a]],validation_genes[[b]]))/
          length(union(discovery_genes[[a]],validation_genes[[b]])))
    }))
  }))
  tab$FDR <- p.adjust(tab$P,"BH")
  write.csv(tab,prior,row.names=FALSE)
  write.csv(data.frame(gene=universe),file.path(out,"NMF_eligible_gene_universe.csv"),row.names=FALSE)
  information <- list(corrected_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
    analysis="Enrichment of frozen discovery signatures within validation de novo NMF programs",
    background_n=length(universe),n_fitted_NMF_models=length(programs),
    all_saved_W_gene_sets_identical=TRUE,
    background_source="rownames(W) from the saved validation NMF fits; not all measured genes",
    measured_gene_count=length(measured_universe),actual_multiNMF_seed=actual_NMF_seed,
    multiple_testing="BH across all 16 discovery-by-validation comparisons",
    jaccard="Original measured-gene raw Jaccard preserved in jaccard; descriptive only",
    additional_overlap="jaccard_eligible uses the same eligible universe as Fisher; unrestricted sets also recorded",
    refitted_NMF=FALSE,
    limitation="Discovery W candidate universe was not reconstructed here; inferential test treats discovery signatures as fixed and conditions on validation-eligible genes. It does not claim a symmetric two-cohort candidate-universe test.",
    prior_full_gene_results="Audit only; their P values must not annotate a figure or support a claim")
  jsonlite::write_json(information,file.path(audit,"NMF_overlap_background_correction.json"),
    pretty=TRUE,auto_unbox=TRUE)
  print(tab)
  invisible(tab)
}

if(sys.nframe()==0L) {
  suppressPackageStartupMessages({library(Seurat);library(jsonlite)})
  args<-commandArgs(trailingOnly=TRUE)
  script_arg<-commandArgs(trailingOnly=FALSE)
  script_file<-sub("^--file=","",script_arg[grepl("^--file=",script_arg)][1])
  root<-if(length(args))args[1] else normalizePath(file.path(dirname(script_file),".."),mustWork=TRUE)
  out<-file.path(root,"data/processed")
  saved<-readRDS(file.path(out,"validation_de_novo_NMF.rds"))
  scored<-readRDS(file.path(out,"validation_primary_scored.rds"))
  discovery<-read.csv(file.path(out,"discovery_NMF_genes.csv"))
  recompute_NMF_overlap(saved$programs,saved$metaprograms$metaprograms.genes,
    split(discovery$gene,discovery$program),rownames(scored),out,file.path(root,"audit"),123L)
}
