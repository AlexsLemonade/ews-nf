#!/bin/bash

# Manually launch a tmux session and run the ewsnf workflow from a local clone.
# Usage: bash launch_ewsnf_testing.sh [/path/to/workflow/repo]

set -euo pipefail

# pin the Nextflow version
NXF_VER=${NXF_VER:-25.10.7}

# directory of the cloned workflow repo
workflow_dir=${1:-${WORKFLOW_DIR:-$HOME/ewsnf}}

# config file location, as a path relative to the repo root
config_file=${CONFIG_FILE:-ccdl_profiles.config}

tmux_session="tmux_ewsnf_testing"

if [[ ! -f "${workflow_dir}/main.nf" ]]; then
  echo "Error: main.nf not found in '${workflow_dir}'." 1>&2
  exit 1
fi

# if the config isn't at the repo root, look for it elsewhere in the repo
if [[ ! -f "${workflow_dir}/${config_file}" ]]; then
  found_config=$(cd "$workflow_dir" && find . -name "$(basename "$config_file")" -not -path './.git/*' -print -quit)
  if [[ -z "$found_config" ]]; then
    echo "Error: '${config_file}' not found in '${workflow_dir}'." 1>&2
    exit 1
  fi
  config_file="${found_config#./}"
fi

echo "Using config: ${config_file}"

if tmux has-session -t "$tmux_session" 2>/dev/null; then
  echo "Error: tmux session '${tmux_session}' already exists. Attach or kill it first." 1>&2
  exit 1
fi

tmux new-session -d -s "$tmux_session" -c "$workflow_dir"

# keep the pane around after nextflow exits so the output is still readable
tmux set-option -t "$tmux_session" remain-on-exit on

tmux set-environment -t "$tmux_session" NXF_VER "$NXF_VER"
tmux set-environment -t "$tmux_session" NXF_JVM_ARGS "-XX:InitialRAMPercentage=10 -XX:MaxRAMPercentage=40"
tmux set-environment -t "$tmux_session" PATH "$HOME/.local/bin:$HOME/bin:$PATH"

# respawn the pane so it picks up the session environment set above
tmux respawn-pane -t "$tmux_session" -k -c "$workflow_dir" \
  "nextflow run main.nf -c ${config_file} -profile scpca_testing -with-tower -resume"

echo "Launched session '${tmux_session}' in ${workflow_dir}."
echo "Attach with: tmux attach -t ${tmux_session}"
