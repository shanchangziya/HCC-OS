#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
analysis_root="$(cd "${script_dir}/.." && pwd)"
repo_root="$(cd "${analysis_root}/../.." && pwd)"
raw_root="${HCC_OS_MODEL_RAW_DIR:-${repo_root}/revision/model_validation/data/raw}"
processed_root="${analysis_root}/data/processed"
logs_root="${analysis_root}/logs"
task_logs="${logs_root}/server_resume_tasks"
max_jobs="${OSARS_MAX_JOBS:-24}"

mkdir -p "${task_logs}"

if find "${raw_root}" -maxdepth 1 -type f -exec basename {} \; | grep -qi 'ICGC'; then
  printf 'Refusing TCGA development run because an ICGC file is present before model lock.\n' >&2
  exit 1
fi

for fold in 1 2 3 4 5; do
  test -s "${processed_root}/TCGA_fold${fold}_input.rds"
done

tasks=()
for fold in 1 2 3 4 5; do
  for config_id in $(seq 1 117); do
    checkpoint="${processed_root}/fold${fold}_checkpoints/config_$(printf '%03d' "${config_id}").rds"
    if [[ ! -s "${checkpoint}" ]]; then
      tasks+=("${fold}:${config_id}")
    fi
  done
done

printf 'Missing fold/configuration tasks at resume: %d\n' "${#tasks[@]}"
printf '%s\n' "${tasks[@]}" > "${logs_root}/server_resume_task_manifest.txt"

failed=0
for task in "${tasks[@]}"; do
  fold="${task%%:*}"
  config_id="${task##*:}"
  while [[ "$(jobs -rp | wc -l | tr -d ' ')" -ge "${max_jobs}" ]]; do
    if ! wait -n; then
      failed=1
    fi
  done
  (
    OSARS_CONFIG_IDS="${config_id}" Rscript "${script_dir}/02_run_tcga_fold.R" "${fold}" \
      > "${task_logs}/fold${fold}_config$(printf '%03d' "${config_id}").console.log" 2>&1
  ) &
done

while [[ "$(jobs -rp | wc -l | tr -d ' ')" -gt 0 ]]; do
  if ! wait -n; then
    failed=1
  fi
done
if [[ "${failed}" -ne 0 ]]; then
  printf 'At least one resumed R process exited non-zero.\n' >&2
  exit 1
fi

checkpoint_count="$(find "${processed_root}" -path '*/fold*_checkpoints/config_*.rds' -type f | wc -l | tr -d ' ')"
if [[ "${checkpoint_count}" -ne 585 ]]; then
  printf 'Expected 585 checkpoints, found %s.\n' "${checkpoint_count}" >&2
  exit 1
fi

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
printf 'Parallel resume, TCGA aggregation, and model lock completed.\n'
