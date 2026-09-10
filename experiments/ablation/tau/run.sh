#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
VENV_DIR="${REPO_ROOT}/.venv"

repo_abs_path() {
  case "$1" in
    /*) printf "%s\n" "$1" ;;
    *) printf "%s/%s\n" "$REPO_ROOT" "$1" ;;
  esac
}

DATA_PATH="${1:-}"
if [ -z "$DATA_PATH" ]; then
    exit 1
fi

OUTPUT_PATH="${2:-}"
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

K="${3:-}"
if [ -z "$K" ]; then
    exit 1
fi

PROPORTIONS_PATH="${4:-}"
if [ -z "$PROPORTIONS_PATH" ]; then
    exit 1
fi

DATA_PATH="$(repo_abs_path "$DATA_PATH")"
OUTPUT_PATH="$(repo_abs_path "$OUTPUT_PATH")"
PROPORTIONS_PATH="$(repo_abs_path "$PROPORTIONS_PATH")"

values=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)

mkdir -p "$OUTPUT_PATH/logs/"

for val in "${values[@]}"
do
(
    cd "$REPO_ROOT/methods/SNMF"
    mkdir -p "$OUTPUT_PATH/v$val/"
    bash run.sh \
        "$DATA_PATH" \
        "$OUTPUT_PATH/v$val/" \
        "$val" \
        "KL_NB" \
        "$K" \
        "$PROPORTIONS_PATH" \
        "42"
  ) > "$OUTPUT_PATH/logs/SNMF_v${val}.log" 2>&1 &

  sleep 10
done

wait

if [ -d "$VENV_DIR" ]; then
  source "$VENV_DIR/bin/activate"
fi
python "$SCRIPT_DIR/plot_metrics.py" \
  "$OUTPUT_PATH" \
  "$PROPORTIONS_PATH"
