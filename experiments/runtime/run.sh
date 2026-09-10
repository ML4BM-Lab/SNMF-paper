#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
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

DATA_PATH="$(repo_abs_path "$DATA_PATH")"
OUTPUT_PATH="$(repo_abs_path "$OUTPUT_PATH")"
if [ -n "$PROPORTIONS_PATH" ]; then
    PROPORTIONS_PATH="$(repo_abs_path "$PROPORTIONS_PATH")"
fi

seeds=(1 2 3 4 5)

mkdir -p "$OUTPUT_PATH/logs"

for seed in "${seeds[@]}"
do
  (
    cd "$SCRIPT_DIR"
    mkdir -p "$OUTPUT_PATH/cpu/$seed/"
    bash cpu.sh \
        "$DATA_PATH" \
        "$OUTPUT_PATH/cpu/$seed/" \
        "0.4" \
        "KL_NB" \
        "$K" \
        "$PROPORTIONS_PATH" \
        "$seed"
  ) > "$OUTPUT_PATH/logs/cpu_${seed}.log" 2>&1 &

  sleep 20

  (
    cd "$REPO_ROOT/methods/SNMF"
    mkdir -p "$OUTPUT_PATH/gpu/$seed/"
    bash run.sh \
        "$DATA_PATH" \
        "$OUTPUT_PATH/gpu/$seed/" \
        "0.4" \
        "KL_NB" \
        "$K" \
        "$PROPORTIONS_PATH" \
        "$seed"
  ) > "$OUTPUT_PATH/logs/gpu_${seed}.log" 2>&1 &

  sleep 20
done

wait

if [ -d "$VENV_DIR" ]; then
  source "$VENV_DIR/bin/activate"
fi
python "$SCRIPT_DIR/plot_time.py" \
  "$OUTPUT_PATH"
