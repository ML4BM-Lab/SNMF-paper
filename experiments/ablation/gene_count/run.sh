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

MARKERS_PATH="${2:-}"
if [ -z "$MARKERS_PATH" ]; then
    exit 1
fi

SNMF_TAU="${3:-}"
if [ -z "$SNMF_TAU" ]; then
    exit 1
fi

OUTPUT_PATH="${4:-}"
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

K="${5:-}"
if [ -z "$K" ]; then
    exit 1
fi

PROPORTIONS_PATH="${6:-}"
if [ -z "$PROPORTIONS_PATH" ]; then
    exit 1
fi

DATA_PATH="$(repo_abs_path "$DATA_PATH")"
MARKERS_PATH="$(repo_abs_path "$MARKERS_PATH")"
OUTPUT_PATH="$(repo_abs_path "$OUTPUT_PATH")"
PROPORTIONS_PATH="$(repo_abs_path "$PROPORTIONS_PATH")"

ngenes=(100 200 400 800 1600 3200 6400 12800 25600)

mkdir -p "$OUTPUT_PATH/tmp"
mkdir -p "$OUTPUT_PATH/logs"
mkdir -p "$SCRIPT_DIR/outputs/logs"

for val in "${ngenes[@]}"
do
  mkdir -p "$OUTPUT_PATH/v$val"

  echo "Subsampling genes..."
  (
    cd "$SCRIPT_DIR"
    sbatch --parsable --wait ./subsample_data.slurm \
      "$DATA_PATH" \
      "$MARKERS_PATH" \
      "$val" \
      "$OUTPUT_PATH"
  )
  echo "Genes subsampled!"

  (
    bash "$REPO_ROOT/experiments/benchmarking/run_benchmark.sh" \
        --data_path="$OUTPUT_PATH/tmp/v$val.csv" \
        --markers_path="$OUTPUT_PATH/tmp/mg$val.csv" \
        --output_path="$OUTPUT_PATH/v$val/" \
        --k="$K" \
        --proportions_path="$PROPORTIONS_PATH" \
        --snmf_tau="$SNMF_TAU" \
        --starfysh_lr=1e-6 \
  ) > "$OUTPUT_PATH/logs/SNMF_v${val}.log" 2>&1 &

  echo "SNMF Launched"

  sleep 100
done

wait

rm -r "$OUTPUT_PATH/tmp"

if [ -d "$VENV_DIR" ]; then
  source "$VENV_DIR/bin/activate"
fi
python "$SCRIPT_DIR/plot_metrics.py" \
  "$OUTPUT_PATH"
