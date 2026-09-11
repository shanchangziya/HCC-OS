#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
analysis_root="$(cd "${script_dir}/.." && pwd)"
repo_root="$(cd "${analysis_root}/../.." && pwd)"
raw_root="${HCC_OS_MODEL_RAW_DIR:-${repo_root}/revision/model_validation/data/raw}"
logs_root="${analysis_root}/logs"
processed_root="${analysis_root}/data/processed"

mkdir -p "${logs_root}" "${processed_root}"

find "${raw_root}" -maxdepth 1 -type f -exec sh -c \
  'for file do printf "%s\t%s bytes\n" "$(basename "$file")" "$(wc -c < "$file" | tr -d " ")"; done' sh {} + \
  | sort > "${logs_root}/server_prelock_raw_input_inventory.txt"
if find "${raw_root}" -maxdepth 1 -type f -exec basename {} \; | grep -qi 'ICGC'; then
  printf 'Refusing TCGA development run because an ICGC file is present before model lock.\n' >&2
  exit 1
fi

Rscript -e 'sessionInfo(); pkgs <- c("survival","randomForestSRC","glmnet","plsRcox","superpc","gbm","CoxBoost","survivalsvm","dplyr","tibble","miscTools","compareC","mixOmics"); cat("\nMODEL PACKAGE VERSIONS\n"); for (p in pkgs) cat(p, as.character(packageVersion(p)), "\n")' \
  > "${logs_root}/server_R_environment.txt" 2>&1

Rscript "${script_dir}/01_prepare_tcga_folds.R" \
  > "${logs_root}/server_01_prepare_tcga_folds.console.log" 2>&1

pids=()
labels=()
for fold in 1 2 3 4 5; do
  chunk_index=0
  for bounds in "1 30" "31 60" "61 90" "91 117"; do
    chunk_index=$((chunk_index + 1))
    read -r first last <<< "${bounds}"
    config_ids="$(seq -s, "${first}" "${last}")"
    log_file="${logs_root}/server_fold${fold}_chunk${chunk_index}.console.log"
    OSARS_CONFIG_IDS="${config_ids}" Rscript "${script_dir}/02_run_tcga_fold.R" "${fold}" \
      > "${log_file}" 2>&1 &
    pids+=("$!")
    labels+=("fold${fold}_chunk${chunk_index}")
  done
done

failed=0
for i in "${!pids[@]}"; do
  if wait "${pids[$i]}"; then
    printf '%s\tOK\n' "${labels[$i]}"
  else
    printf '%s\tFAILED\n' "${labels[$i]}" >&2
    failed=1
  fi
done
if [[ "${failed}" -ne 0 ]]; then
  exit 1
fi

# Re-open the complete checkpoint set once per fold to write complete fold summaries.
for fold in 1 2 3 4 5; do
  Rscript "${script_dir}/02_run_tcga_fold.R" "${fold}" \
    > "${logs_root}/server_fold${fold}_checkpoint_audit.console.log" 2>&1
done

Rscript "${script_dir}/03_aggregate_select_and_lock.R" \
  > "${logs_root}/server_03_aggregate_select.console.log" 2>&1

Rscript "${script_dir}/04_fit_and_lock_final_tcga_model.R" \
  > "${logs_root}/server_04_lock_final_TCGA_model.console.log" 2>&1

test -s "${processed_root}/MODEL_LOCKED.sha256"
date -u '+%Y-%m-%dT%H:%M:%SZ' > "${processed_root}/SERVER_TCGA_DEVELOPMENT_COMPLETE"

printf 'TCGA-only model development and locking completed. ICGC was not available in this workspace.\n'
