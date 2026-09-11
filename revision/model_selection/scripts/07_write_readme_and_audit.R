#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
script_path <- normalizePath(script_arg[[1]], mustWork = TRUE)
analysis_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(analysis_root, "scripts", "common.R"))

performance <- read.csv(file.path(analysis_root, "all_117_models_5fold_performance.csv"), check.names = FALSE)
icgc <- read.csv(file.path(analysis_root, "ICGC_retrospective_evaluation.csv"), check.names = FALSE)
feature_folds <- read.csv(file.path(analysis_root, "feature_selection_by_fold.csv"), check.names = FALSE)
final_features <- read.csv(file.path(processed_root, "final_full_TCGA_feature_selection.csv"), check.names = FALSE)
locked_bundle <- readRDS(file.path(processed_root, "final_OSARS_model_locked.rds"))
duplicates <- read.csv(file.path(processed_root, "duplicate_predictor_audit.csv"), check.names = FALSE)
historical_duplicates <- read.csv(file.path(processed_root, "historical_duplicate_predictor_audit.csv"), check.names = FALSE)
selected <- performance[performance$selected_final_configuration, , drop = FALSE]
near <- performance[performance$within_0_01_of_best, , drop = FALSE]
if (nrow(selected) != 1L) stop("Expected exactly one selected configuration")

tcga_scores <- read.csv(file.path(processed_root, "final_model_TCGA_apparent_scores.csv"), check.names = FALSE)
tcga_apparent <- harrell_c(tcga_scores$OS.time, tcga_scores$OS, tcga_scores$locked_risk_score)
fold_values <- as.numeric(selected[1, paste0("fold", 1:5, "_Cindex")])
fold_range <- range(fold_values)
stability_text <- if (selected$SD < 0.03) {
  "fold间波动较小"
} else if (selected$SD < 0.05) {
  "存在中等程度的fold间波动"
} else {
  "存在较明显的fold间波动"
}

same_as_historical <- selected$configuration %in% c("StepCox[forward] + GBM", "GBM")
all_complete <- all(performance$five_fold_predictions_complete)
old_tcga <- 0.870493991989319
old_icgc <- 0.766423357664234

prelock_inventory_file <- file.path(logs_root, "server_prelock_raw_input_inventory.txt")
postlock_transfer_file <- file.path(logs_root, "server_postlock_ICGC_transfer_manifest.txt")
prelock_inventory <- if (file.exists(prelock_inventory_file)) {
  paste(readLines(prelock_inventory_file, warn = FALSE), collapse = "; ")
} else {
  "not available"
}
postlock_transfer <- if (file.exists(postlock_transfer_file)) {
  x <- read.delim(postlock_transfer_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  stats::setNames(x[[2]], x[[1]])
} else {
  character()
}
lock_time_utc <- if ("TCGA_development_completed_utc" %in% names(postlock_transfer)) {
  postlock_transfer[["TCGA_development_completed_utc"]]
} else {
  "not available"
}
icgc_transfer_time_utc <- if ("ICGC_transfer_utc" %in% names(postlock_transfer)) {
  postlock_transfer[["ICGC_transfer_utc"]]
} else {
  "not available"
}

near_table <- c(
  "| Rank | Configuration | Mean C-index | SD | Mean input features | Duplicate group |",
  "|---:|---|---:|---:|---:|---|",
  vapply(seq_len(nrow(near)), function(i) sprintf(
    "| %s | %s | %.4f | %.4f | %.1f | %s |",
    near$rank[i], near$configuration[i], near$mean_Cindex[i], near$SD[i],
    near$mean_n_features[i], ifelse(is.na(near$exact_duplicate_group[i]), "", near$exact_duplicate_group[i])
  ), character(1))
)

readme <- c(
  "# OSARS 117-configuration TCGA-only model-selection audit",
  "",
  "## Scope and headline result",
  "",
  paste0(
    "The frozen modelling grid contains **117 configurations (20 single algorithms + 97 two-stage configurations)**, ",
    "not 101. Model development and algorithm selection were repeated using stratified 5-fold cross-validation entirely within TCGA-LIHC. ",
    "The selected configuration was **", selected$configuration, "**, with a mean held-out-fold Harrell C-index of **",
    sprintf("%.4f", selected$mean_Cindex), "** (SD ", sprintf("%.4f", selected$SD), ")."
  ),
  "",
  "ICGC-LIRI-JP was opened only after the TCGA-selected model had been refitted on complete TCGA and saved with a SHA256 lock. Because the historical analysis had already inspected ICGC while choosing the original model, its present result is labelled **retrospective external evaluation**, not untouched or fully independent validation.",
  "",
  "## Frozen inputs and cohort definition",
  "",
  "- TCGA expression: `TCGA_unclipped_448_expression.csv` (343 patients; 448 frozen candidate genes).",
  "- TCGA outcome: `TCGA-LIHC_frozen_risks.csv` (343 patients; 124 deaths/events). Only `ID`, `OS.time`, and `OS` were used; the historical risk score was not used.",
  "- ICGC retrospective evaluation: `ICGC_original_448_expression.csv` (243 patients).",
  "- Patient inclusion, survival-time/status definitions, and expression matrices were not redefined.",
  "- Input file SHA256 hashes are recorded in `data/processed/TCGA_input_manifest.csv` and the external-evaluation log.",
  paste0("- Server pre-lock raw-input inventory: ", prelock_inventory, "."),
  paste0(
    "- The TCGA development run was locked at ", lock_time_utc,
    "; the ICGC file was transferred into the isolated server workspace only afterwards, at ",
    icgc_transfer_time_utc, "."
  ),
  "",
  "## Revised leakage-controlled workflow",
  "",
  "1. A fixed master seed (`5201314`) generated event-stratified 5-fold assignments. Held-out folds contain 67–69 patients and event fractions of 35.8%–36.2%.",
  "2. For each outer fold, 1st/99th percentile winsorization bounds, mean-imputation values, and near-zero-variance filtering were estimated using only the 80% training partition and then applied unchanged to its held-out partition.",
  "3. Univariable Cox screening (`P < 0.01`) was repeated independently inside each training partition. The five folds retained ",
  paste0("   ", paste(feature_folds$genes_after_univariable_cox_p_lt_0_01, collapse = ", "), " genes, respectively."),
  "4. All model fitting and all algorithm-specific tuning (10-fold internal tuning where used by the frozen Mime implementation) used only that outer training partition.",
  "5. Risk-score orientation was determined from training-partition predictions only and then fixed before scoring the held-out patients. Held-out outcomes were used only to calculate Harrell's C-index.",
  "6. The 117 configurations were ranked by mean outer-fold C-index. ICGC was not used as a tie-breaker.",
  "7. The pre-specified tie rule was: highest mean C-index; exact numerical ties then favour fewer input features, then a single-stage configuration, then frozen configuration order.",
  "8. The chosen configuration was refitted once on complete TCGA, serialized, and SHA256-locked. Only then was the external-evaluation script allowed to open ICGC.",
  "",
  "## Leakage audit",
  "",
  "| Potential leakage | Revised analysis |",
  "|---|---|",
  "| Univariable Cox on all TCGA before CV | Prevented; repeated in each outer training fold |",
  "| Global TCGA winsorization/imputation before CV | Prevented; parameters estimated in each outer training fold |",
  "| Held-out fold used for tuning | Prevented; tuning calls receive training-fold data only |",
  "| Risk direction inferred from held-out outcomes | Prevented; orientation fixed from training predictions |",
  "| Model changed after viewing held-out performance | Prevented; one frozen grid and one pre-specified ranking rule |",
  "| ICGC used for feature/algorithm/tuning/cutoff decisions | Prevented in revised scripts; model hash verified before ICGC is opened |",
  "",
  "No external cohort path or object is referenced by the fold-preparation, fold-fitting, aggregation/selection, or complete-TCGA locking scripts.",
  "",
  "## Original versus revised workflow",
  "",
  "| Item | Historical workflow | Revised workflow |",
  "|---|---|---|",
  "| Configuration count | Manuscript stated 101; frozen object contains 117 | Audited and reports all 117 |",
  "| Univariable screen | Complete TCGA before the reported multi-model analysis | Repeated inside each outer training fold |",
  "| Expression capping/imputation | Estimated cohort-wise before model evaluation | Estimated only in each TCGA training fold; full-TCGA values locked for later use |",
  "| Algorithm selection | ICGC results were available and used when prioritising the historical StepCox[forward] + GBM/GBM predictor | Mean 5-fold CV C-index within TCGA only |",
  "| External cohort label | Previously treated as independent validation | Retrospective external evaluation |",
  "| Risk direction | Historical evaluator fitted `Surv ~ RS` separately in each evaluated cohort | Fixed using training data before held-out scoring |",
  "",
  "## Duplicate-predictor audit",
  "",
  paste0(
    "The frozen historical object contains ", nrow(historical_duplicates),
    " exact or rank-equivalent predictor pairs. In particular, `StepCox[forward] + GBM` and `GBM` were exactly identical in the historical TCGA and ICGC score vectors. ",
    "This is expected from the implementation because StepCox `forward` starts from the full Cox model and therefore may retain the complete input set."
  ),
  paste0(
    "In the revised out-of-fold predictions, ", nrow(duplicates),
    " exact or rank-equivalent pairs were detected; they are retained as separate configurations and marked in `all_117_models_5fold_performance.csv`."
  ),
  "",
  "## Models within 0.01 of the best",
  "",
  near_table,
  "",
  "## Final locked model",
  "",
  paste0("- Configuration: `", selected$configuration, "`."),
  paste0("- Mean 5-fold C-index within TCGA: ", sprintf("%.4f", selected$mean_Cindex), "."),
  paste0("- Fold values: ", paste(sprintf("%.4f", fold_values), collapse = ", "), "."),
  paste0("- Fold SD: ", sprintf("%.4f", selected$SD), "; range ", sprintf("%.4f", fold_range[1]), "–", sprintf("%.4f", fold_range[2]), " (", stability_text, ")."),
  paste0("- Complete-TCGA univariable genes: ", final_features$genes_after_univariable_cox_p_lt_0_01, "."),
  paste0("- Final model input features: ", final_features$final_model_input_features, "."),
  "- Reproducible prediction is implemented in `scripts/predict_locked_model.R`; preprocessing vectors, feature order, model object, parameters, and risk orientation are all stored in `data/processed/final_OSARS_model_locked.rds`.",
  "- The explicit nine-component PLS-Cox calculation and its parameter tables are exported in `final_model_risk_formula.txt` and `data/processed/final_model_*`. These reproduce serialized predictions to machine precision without changing the model hash.",
  "",
  "## ICGC retrospective external evaluation",
  "",
  paste0(
    "The locked model yielded a Harrell C-index of **", sprintf("%.4f", icgc$Cindex),
    "** (95% CI ", sprintf("%.4f", icgc$CI95_low), "–", sprintf("%.4f", icgc$CI95_high),
    "; n=", icgc$n_patients, ", events=", icgc$n_events, ")."
  ),
  "No ICGC feature selection, algorithm selection, hyperparameter tuning, risk-direction selection, cutoff optimisation, or model modification was performed.",
  "",
  "## Interpretation of performance change",
  "",
  paste0(
    "The historical TCGA value (0.8705) was an apparent in-sample C-index, whereas the revised primary TCGA estimate is held-out-fold CV performance (",
    sprintf("%.4f", selected$mean_Cindex), "). Their difference (", sprintf("%+.4f", selected$mean_Cindex - old_tcga),
    ") is not a like-for-like degradation; it mainly reflects removal of resubstitution optimism and fold-wise repetition of feature selection."
  ),
  paste0(
    "The historical ICGC C-index was 0.7664 and the locked revised model gives ", sprintf("%.4f", icgc$Cindex),
    " (difference ", sprintf("%+.4f", icgc$Cindex - old_icgc), "). This is a retrospective comparison because ICGC had been examined historically."
  ),
  paste0("For reference only, the refitted locked model's complete-TCGA apparent C-index is ", sprintf("%.4f", tcga_apparent), "."),
  "",
  "## Statistical limitation",
  "",
  "The requested five outer folds provide honest held-out predictions for every individual configuration and prevent fold-level feature/tuning leakage. However, the same set of five fold estimates is used to rank 117 configurations and to quote the selected configuration's mean; therefore some winner's-curse/model-selection optimism can remain. These values are explicitly called **internal cross-validation/model-selection performance**, not independent test performance. A fully nested outer assessment of the entire 117-configuration selection strategy would be needed for an unbiased post-selection performance estimate.",
  "",
  "## Software fidelity and recorded deviations",
  "",
  "- The modelling code vendors the upstream `ML.Dev.Prog.Sig` implementation used to reconstruct the frozen grid and keeps the original algorithm parameters (RSF 1,000 trees/nodesize 10; internal 10-fold tuning; GBM up to 10,000 trees; CoxBoost penalty/step tuning; and the published StepCox directions).",
  "- Unused plotting-only package imports were removed from the vendored function; this does not alter model calculations.",
  "- SuperPC's upstream infinite retry loop was bounded to three attempts so a failed configuration is reported rather than hanging indefinitely.",
  paste0(
    "- Server-side package versions were recorded in the locked object: ",
    paste(names(locked_bundle$package_versions), locked_bundle$package_versions, sep = " ", collapse = "; "), "."
  ),
  "- GBM was run with one worker to keep parallel execution deterministic; its statistical parameters are unchanged.",
  "",
  "## Output map",
  "",
  "- `all_117_models_5fold_performance.csv`: ranked fold-level and summary C-indices, CI, completion, feature counts, and duplicate groups.",
  "- `cv_patient_predictions.csv`: long-format out-of-fold risk scores for every patient/configuration.",
  "- `feature_selection_by_fold.csv`: train-only preprocessing/univariable screening counts.",
  "- `final_model_summary.txt`: locked algorithm, CV performance, parameters, and risk-score method.",
  "- `final_model_risk_formula.txt` and `data/processed/final_model_*`: explicit preprocessing, PLS projection, Cox coefficients, and formula validation.",
  "- `ICGC_retrospective_evaluation.csv`: locked-model retrospective external result.",
  "- `data/processed/duplicate_predictor_audit.csv`: revised exact/rank-equivalent pairs.",
  "- `figures/`: 17-cm-wide combined and individual PDF/SVG/PNG panels.",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
)
writeLines(readme, file.path(analysis_root, "README.md"))

answers <- c(
  paste0("1. 117 configurations all completed five-fold CV: ", ifelse(all_complete, "Yes", paste0("No; ", sum(performance$five_fold_predictions_complete), "/117 were complete")), "."),
  paste0("2. Rank 1 in TCGA CV: ", selected$configuration, "."),
  paste0("3. Is it the historical StepCox[forward] + GBM / GBM predictor family: ", ifelse(same_as_historical, "Yes", "No"), "."),
  paste0("4. Mean five-fold C-index: ", sprintf("%.6f", selected$mean_Cindex), "."),
  paste0("5. Fold stability: SD=", sprintf("%.6f", selected$SD), ", range=", sprintf("%.6f", fold_range[1]), "–", sprintf("%.6f", fold_range[2]), "; ", stability_text, "."),
  paste0("6. Models with delta <0.01: ", nrow(near), "."),
  paste0("7. Clear superiority versus near-equivalent alternatives: ", ifelse(nrow(near) == 1L, "Only the best lies within 0.01", paste0(nrow(near), " models lie within 0.01; multiple near-equivalent alternatives exist")), "."),
  paste0("8. Locked-model ICGC retrospective C-index: ", sprintf("%.6f", icgc$Cindex), " (95% CI ", sprintf("%.6f", icgc$CI95_low), "–", sprintf("%.6f", icgc$CI95_high), ")."),
  paste0("9. Change versus historical values: TCGA primary estimate changed from apparent 0.870494 to CV ", sprintf("%.6f", selected$mean_Cindex), " (not like-for-like); ICGC changed from 0.766423 to ", sprintf("%.6f", icgc$Cindex), ".")
)
writeLines(answers, file.path(analysis_root, "audit_conclusions.txt"))
cat(paste(answers, collapse = "\n"), "\n")
