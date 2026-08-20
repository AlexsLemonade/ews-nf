#!/bin/bash

# Runs the ews-nf workflow. Intended to be launched inside tmux by
# launch_ewsnf_testing.sh, which sets WORKFLOW_DIR / CONFIG_FILE / NXF_VER.

set -uo pipefail

workflow_dir=${WORKFLOW_DIR:?WORKFLOW_DIR not set}
config_file=${CONFIG_FILE:-ccdl_profiles.config}
datetime=$(date "+%Y-%m-%dT%H%M")

cd "$workflow_dir"

run_name="ewsnf_testing_${datetime}"
log_dir="${workflow_dir}/run_logs/${datetime}"
mkdir -p "$log_dir"

#### Tower credentials ####
# -with-tower silently does nothing useful without a token, so fetch it the same
# way run_nextflow.sh does. Set SKIP_TOWER=true to run without it.
if [[ "${SKIP_TOWER:-false}" != "true" ]]; then
  aws_secrets=$(aws secretsmanager get-secret-value \
    --secret-id scpca_service_account_token | jq -r '.SecretString') || aws_secrets=""

  if [[ -n "$aws_secrets" ]]; then
    OP_SERVICE_ACCOUNT_TOKEN=$(jq -r '.op_token' <<< "$aws_secrets")
    export OP_SERVICE_ACCOUNT_TOKEN
    TOWER_ACCESS_TOKEN=$(op read "$(jq -r '.op_seqera_token' <<< "$aws_secrets")")
    export TOWER_ACCESS_TOKEN
  fi

  if [[ -z "${TOWER_ACCESS_TOKEN:-}" ]]; then
    echo "WARNING: TOWER_ACCESS_TOKEN is empty; Tower reporting will not work." 1>&2
  else
    echo "Tower token loaded (${#TOWER_ACCESS_TOKEN} chars)."
  fi
fi

#### Run the workflow ####

echo "Run name: ${run_name}"
echo "Logs:     ${log_dir}"
echo ""

# -ansi-log false keeps the tee'd log readable; trace/report give per-task detail
nextflow run main.nf \
  -c "$config_file" \
  -profile scpca_testing \
  -name "$run_name" \
  -ansi-log false \
  -with-report "${log_dir}/report.html" \
  -with-trace "${log_dir}/trace.txt" \
  -with-tower \
  ${NEXTFLOW_OPTIONS:-} \
  2>&1 | tee "${log_dir}/nextflow_stdout.log"

exit_code=${PIPESTATUS[0]}

# keep the full .nextflow.log alongside the rest, masking the token if present
if [[ -f .nextflow.log ]]; then
  cp .nextflow.log "${log_dir}/nextflow.log"
  if [[ -n "${TOWER_ACCESS_TOKEN:-}" ]]; then
    sed -i "s/${TOWER_ACCESS_TOKEN}/<TOWER_ACCESS_TOKEN>/g" "${log_dir}"/*.log
  fi
fi

#### Collect logs from failed tasks ####

if [[ $exit_code -ne 0 ]]; then
  echo ""
  echo "=============================================="
  echo "Run failed (exit ${exit_code}). Failed tasks:"
  echo "=============================================="

  nextflow log "$run_name" -f process,tag,exit,workdir -F "status == 'FAILED'" \
    | tee "${log_dir}/failed_tasks.txt"

  # pull the task logs down from the S3 work dir so they can be read locally
  while read -r process tag task_exit workdir; do
    [[ -z "${workdir:-}" ]] && continue
    dest="${log_dir}/failed/${process}"
    mkdir -p "$dest"
    echo "Fetching logs for ${process} (${tag}) from ${workdir}"
    for f in .command.err .command.out .command.log .command.sh .exitcode; do
      if [[ "$workdir" == s3://* ]]; then
        aws s3 cp "${workdir}/${f}" "${dest}/${f}" --quiet 2>/dev/null
      else
        cp "${workdir}/${f}" "${dest}/${f}" 2>/dev/null
      fi
    done
  done < "${log_dir}/failed_tasks.txt"

  echo ""
  echo "Task logs saved under ${log_dir}/failed/"
  echo "Start with: cat ${log_dir}/failed/*/.command.err"
fi

echo ""
echo "Finished with exit code ${exit_code}. Full logs in ${log_dir}"
echo "Pane will stay open; press Ctrl-b d to detach."
