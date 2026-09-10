#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

repo_abs_path() {
  case "$1" in
    /*) printf "%s\n" "$1" ;;
    *) printf "%s/%s\n" "$REPO_ROOT" "$1" ;;
  esac
}

echo "[4] Starfysh"

DATA_PATH="${1:-}"
if [ -z "$DATA_PATH" ]; then
    exit 1
fi

MARKERS_PATH="${2:-}"
if [ -z "$MARKERS_PATH" ]; then
    exit 1
fi

OUTPUT_PATH="${3:-}"
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

LEARNING_RATE="${4:-}"
if [ -z "$LEARNING_RATE" ]; then
    exit 1
fi

SEED="${5:-}"
if [ -z "$SEED" ]; then
    exit 1
fi

DATA_PATH="$(repo_abs_path "$DATA_PATH")"
MARKERS_PATH="$(repo_abs_path "$MARKERS_PATH")"
OUTPUT_PATH="$(repo_abs_path "$OUTPUT_PATH")"

cd "$SCRIPT_DIR"
mkdir -p "$SCRIPT_DIR/logs"

# Make temporary folder
mkdir -p "$OUTPUT_PATH/tmp"

echo "Loading data..."
jid1=$(sbatch --parsable --wait ./load_data.slurm $DATA_PATH $MARKERS_PATH $OUTPUT_PATH)
echo "Data loaded!"

echo "Deconvolution started..."
jid2=$(sbatch --parsable --wait ./starfysh.slurm $OUTPUT_PATH $LEARNING_RATE $SEED)
echo "Deconvolution finished!"

sacct -j $jid1 --format=JobID,JobName,MaxRSS,Elapsed,State > $OUTPUT_PATH/preprocessing_sacct.log
sacct -j $jid2 --format=JobID,JobName,MaxRSS,Elapsed,State > $OUTPUT_PATH/sacct.log

rm -r $OUTPUT_PATH/tmp

echo "[4] Starfysh finished"
echo
