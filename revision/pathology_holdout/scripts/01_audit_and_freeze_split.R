#!/usr/bin/env Rscript

# Reviewer 2 Major Comment 6: source audit and one-time strict patient split.
# This script deliberately stops after freezing the split. It does not fit or
# evaluate a pathology model.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 3:5) {
  stop("Usage: Rscript 01_audit_and_freeze_split.R SOURCE_DIR CLINICAL_CSV FRESH_OUTPUT_DIR [ORIGINAL_FEATURE_SCRIPT] [ORIGINAL_MODEL_SCRIPT]")
}

source_dir <- normalizePath(args[[1]], mustWork = TRUE)
clinical_file <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- path.expand(args[[3]])
provenance_scripts <- if(length(args)>3)normalizePath(args[4:length(args)],mustWork=TRUE) else character()

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
split_file <- file.path(output_dir, "pathology_strict_holdout_split.csv")
if (file.exists(split_file)) {
  stop("Frozen split already exists; refusing to overwrite or regenerate: ", split_file)
}

required_sources <- c(
  file.path(source_dir, "Resnet.Rdata"),
  file.path(source_dir, "pathdat.Rdata"),
  file.path(source_dir, "resnet50_features.csv"),
  clinical_file,
  provenance_scripts
)
if (any(!file.exists(required_sources))) {
  stop("Missing source(s): ", paste(required_sources[!file.exists(required_sources)], collapse = "; "))
}

extract_patient_id <- function(x) {
  hit <- regexpr("TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  ok <- hit > 0L
  out[ok] <- regmatches(x, hit)[ok]
  out
}

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "")
}

# Load the already extracted patient-level features and frozen molecular score.
resnet_env <- new.env(parent = emptyenv())
load(file.path(source_dir, "Resnet.Rdata"), envir = resnet_env)
stopifnot(all(c("resnet2", "rp", "resnet_RS_sig") %in% ls(resnet_env)))
patient_features <- resnet_env$resnet2
frozen_molecular <- resnet_env$rp
expected_features <- paste0("resnet", 0:2047)
stopifnot(
  is.data.frame(patient_features),
  identical(colnames(patient_features), expected_features),
  ncol(patient_features) == 2048L,
  !anyDuplicated(rownames(patient_features))
)
feature_matrix_all <- as.matrix(patient_features)
storage.mode(feature_matrix_all) <- "double"
stopifnot(all(is.finite(feature_matrix_all)))

# Re-audit tile and WSI provenance without re-extracting any ResNet features.
tile_table <- read.csv(
  file.path(source_dir, "resnet50_features.csv"),
  header = TRUE, row.names = 1, check.names = FALSE
)
stopifnot(ncol(tile_table) == 2048L)
colnames(tile_table) <- expected_features
tile_matrix <- as.matrix(tile_table)
storage.mode(tile_matrix) <- "double"
stopifnot(all(is.finite(tile_matrix)))

tile_paths <- rownames(tile_table)
tile_patient <- extract_patient_id(tile_paths)
tile_slide <- sub("^.*/([^/]+\\.svs)/score/.*$", "\\1", tile_paths, perl = TRUE)
stopifnot(!anyNA(tile_patient), all(grepl("\\.svs$", tile_slide)))

tile_counts <- table(tile_patient)
patient_sum <- rowsum(tile_matrix, group = tile_patient, reorder = FALSE)
patient_mean_from_tiles <- patient_sum / as.numeric(tile_counts[rownames(patient_sum)])
saved_order <- rownames(patient_features)
stopifnot(setequal(rownames(patient_mean_from_tiles), saved_order))
aggregation_max_abs_difference <- max(abs(
  patient_mean_from_tiles[saved_order, expected_features, drop = FALSE] -
    feature_matrix_all[saved_order, expected_features, drop = FALSE]
))
stopifnot(aggregation_max_abs_difference < 1e-6)

patient_slide_pairs <- unique(data.frame(
  patient_id = tile_patient,
  slide_id = tile_slide,
  stringsAsFactors = FALSE
))
patient_tile_audit <- data.frame(
  patient_id = names(tile_counts),
  n_tiles = as.integer(tile_counts),
  n_wsi = as.integer(table(patient_slide_pairs$patient_id)[names(tile_counts)]),
  stringsAsFactors = FALSE
)
patient_tile_audit <- patient_tile_audit[order(patient_tile_audit$patient_id), ]
write_csv(patient_tile_audit, file.path(output_dir, "pathology_patient_tile_counts.csv"))

# Collapse molecular aliquot IDs only after confirming within-patient agreement.
molecular_patient <- extract_patient_id(as.character(frozen_molecular$ID))
stopifnot(!anyNA(molecular_patient))
molecular_raw <- data.frame(
  patient_id = molecular_patient,
  sample_id = as.character(frozen_molecular$ID),
  survival_time_days = as.numeric(frozen_molecular$OS.time),
  event = as.integer(frozen_molecular$OS),
  OSARS_score = as.numeric(frozen_molecular$RS),
  stringsAsFactors = FALSE
)

duplicate_counts <- table(molecular_raw$patient_id)
duplicate_ids <- names(duplicate_counts)[duplicate_counts > 1L]
duplicate_audit <- data.frame(
  patient_id = duplicate_ids,
  n_molecular_rows = as.integer(duplicate_counts[duplicate_ids]),
  stringsAsFactors = FALSE
)
write_csv(duplicate_audit, file.path(output_dir, "pathology_molecular_duplicate_patients.csv"))

if (length(duplicate_ids)) {
  for (id in duplicate_ids) {
    d <- molecular_raw[molecular_raw$patient_id == id, , drop = FALSE]
    for (variable in c("survival_time_days", "event", "OSARS_score")) {
      if (length(unique(d[[variable]][is.finite(d[[variable]])])) != 1L) {
        stop("Discordant molecular duplicates for ", id, " in ", variable)
      }
    }
  }
}
molecular <- molecular_raw[!duplicated(molecular_raw$patient_id), , drop = FALSE]

# Eligibility is determined before splitting and uses no pathology-model outcome.
eligible <- molecular[
  molecular$patient_id %in% rownames(patient_features) &
    is.finite(molecular$survival_time_days) & molecular$survival_time_days > 0 &
    molecular$event %in% c(0L, 1L) & is.finite(molecular$OSARS_score),
  , drop = FALSE
]
eligible <- eligible[order(eligible$patient_id), ]
stopifnot(!anyDuplicated(eligible$patient_id))

clinical <- read.csv(clinical_file, check.names = FALSE, stringsAsFactors = FALSE)
clinical$patient_id <- extract_patient_id(as.character(clinical$ID))
stopifnot(!anyNA(clinical$patient_id), !anyDuplicated(clinical$patient_id))
clinical_keep <- intersect(
  c("patient_id", "Stage", "Stage_Group", "Grade", "Grade_Group", "Sex"),
  names(clinical)
)
eligible <- merge(
  eligible, clinical[, clinical_keep, drop = FALSE],
  by = "patient_id", all.x = TRUE, sort = FALSE
)
eligible <- eligible[match(sort(eligible$patient_id), eligible$patient_id), ]
stopifnot(!anyDuplicated(eligible$patient_id))

eligible_features <- feature_matrix_all[eligible$patient_id, expected_features, drop = FALSE]
stopifnot(nrow(eligible_features) == nrow(eligible), all(is.finite(eligible_features)))

analysis_input <- list(
  patient_data = eligible,
  feature_matrix = eligible_features,
  feature_names = expected_features,
  aggregation = "Arithmetic mean of every available tile from every WSI belonging to a patient",
  screening_target = "Frozen continuous OSARS risk score (RS)",
  original_screening = list(method = "Pearson", abs_r_cutoff = 0.2, nominal_p_cutoff = 0.05),
  original_lasso = list(
    framework = "cv.glmnet Cox LASSO", alpha = 1, nfolds = 10,
    standardize = TRUE, lambda_rule = "lambda.min"
  )
)
saveRDS(analysis_input, file.path(output_dir, "pathology_analysis_input.rds"))

excluded_feature_patients <- setdiff(rownames(patient_features), eligible$patient_id)
excluded_molecular_patients <- setdiff(molecular$patient_id, eligible$patient_id)
excluded <- rbind(
  data.frame(
    patient_id = excluded_feature_patients,
    source = "patient-level ResNet50 features",
    reason = "No complete frozen OSARS plus valid OS time/event match",
    stringsAsFactors = FALSE
  ),
  data.frame(
    patient_id = excluded_molecular_patients,
    source = "frozen molecular cohort",
    reason = "No patient-level ResNet50 feature match",
    stringsAsFactors = FALSE
  )
)
excluded <- unique(excluded)
write_csv(excluded, file.path(output_dir, "pathology_excluded_patients.csv"))

source_audit <- data.frame(
  metric = c(
    "raw_tile_count", "unique_WSI_count", "unique_pathology_patient_count",
    "patient_level_feature_rows", "patient_level_feature_columns",
    "patient_level_feature_duplicate_ID_count", "duplicate_tile_path_count",
    "frozen_molecular_patient_count", "eligible_matched_patient_count",
    "eligible_event_count", "molecular_duplicate_patient_count",
    "patients_with_multiple_WSI", "tile_to_patient_mean_max_abs_difference",
    "clinical_stage_available_count"
  ),
  value = c(
    nrow(tile_table), length(unique(tile_slide)), length(unique(tile_patient)),
    nrow(patient_features), ncol(patient_features), sum(duplicated(rownames(patient_features))),
    sum(duplicated(tile_paths)), nrow(molecular), nrow(eligible),
    sum(eligible$event), length(duplicate_ids), sum(patient_tile_audit$n_wsi > 1L),
    aggregation_max_abs_difference,
    sum(!is.na(eligible$Stage_Group) & nzchar(eligible$Stage_Group))
  ),
  detail = c(
    "Existing ResNet50 tile-feature rows; no feature extraction rerun",
    "Unique .svs identifiers parsed from tile paths",
    "Unique TCGA patient IDs parsed from tile paths",
    "Existing resnet2 object",
    "resnet0 through resnet2047",
    "Duplicate row names in the archived patient-level feature matrix",
    "Duplicate full tile paths in the archived tile-feature matrix",
    "Unique patient IDs in frozen rp object",
    "Complete patient-level features plus frozen OSARS plus valid OS",
    "Deaths among eligible patients",
    "All duplicates would require concordant OS/OSARS values before collapse",
    "Patients contributing tiles from more than one WSI",
    "Fresh arithmetic tile means versus archived resnet2",
    "Stage was not required for eligibility"
  ),
  stringsAsFactors = FALSE
)
write_csv(source_audit, file.path(output_dir, "pathology_source_audit.csv"))

source_manifest <- data.frame(
  source = basename(required_sources),
  bytes = file.info(required_sources)$size,
  md5 = unname(tools::md5sum(required_sources)),
  stringsAsFactors = FALSE
)
write_csv(source_manifest, file.path(output_dir, "pathology_source_manifest.csv"))

# One and only one prespecified split. Sorting IDs above makes the result stable
# to source-row order. Allocation is stratified once by OS event status.
set.seed(20260906)
train_ids <- unlist(lapply(c(0L, 1L), function(event_value) {
  ids <- eligible$patient_id[eligible$event == event_value]
  sample(ids, size = round(0.70 * length(ids)), replace = FALSE)
}), use.names = FALSE)

split_table <- eligible[, c("patient_id", "survival_time_days", "event"), drop = FALSE]
split_table$split <- ifelse(split_table$patient_id %in% train_ids, "train", "test")
split_table <- split_table[, c("patient_id", "split", "survival_time_days", "event")]
split_table <- split_table[order(split_table$patient_id), ]

stopifnot(
  !anyDuplicated(split_table$patient_id),
  length(intersect(
    split_table$patient_id[split_table$split == "train"],
    split_table$patient_id[split_table$split == "test"]
  )) == 0L,
  nrow(split_table) == nrow(eligible),
  all(table(split_table$event, split_table$split) > 0L)
)

# Atomic first write; the file-existence guard above prohibits regeneration.
temporary_split <- tempfile(pattern = "strict_split_", tmpdir = output_dir, fileext = ".csv")
write_csv(split_table, temporary_split)
if (!file.rename(temporary_split, split_file)) {
  stop("Could not atomically freeze split")
}

split_provenance <- data.frame(
  seed = 20260906L,
  split_level = "patient",
  allocation = "single 70/30 event-stratified split",
  patient_order_before_sampling = "lexicographic patient_id",
  training_n = sum(split_table$split == "train"),
  testing_n = sum(split_table$split == "test"),
  training_events = sum(split_table$event[split_table$split == "train"]),
  testing_events = sum(split_table$event[split_table$split == "test"]),
  train_test_overlap = 0L,
  split_md5 = unname(tools::md5sum(split_file)),
  seed_retry = "No",
  stringsAsFactors = FALSE
)
write_csv(split_provenance, file.path(output_dir, "pathology_split_provenance.csv"))

cat("FROZEN SPLIT CREATED; MODELING NOT STARTED\n")
print(split_provenance)
