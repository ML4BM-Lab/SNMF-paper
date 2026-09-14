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

echo "[10] FAST:"

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

SEED="${4:-}"
if [ -z "$SEED" ]; then
    exit 1
fi

PROPORTIONS_PATH="${5:-}"
if [ -z "$PROPORTIONS_PATH" ]; then
    exit 1
fi

DATA_PATH="$(repo_abs_path "$DATA_PATH")"
OUTPUT_PATH="$(repo_abs_path "$OUTPUT_PATH")"

cd "$SCRIPT_DIR"
mkdir -p "$SCRIPT_DIR/logs"

# Load R
module purge
module load R/4.4.1-gfbf-2023a

# Make temporary folder
mkdir -p "$OUTPUT_PATH/tmp"

echo "Loading data..."
jid1=$(sbatch --parsable --wait ./load_data.slurm $DATA_PATH $OUTPUT_PATH $K)
echo "Data loaded!"

echo "Deconvolution started..."
jid2=$(sbatch --parsable --wait ./FAST.slurm $OUTPUT_PATH $SEED)
echo "Deconvolution finished!"

echo "Annotating cell types..."
Rscript ./annotate.R $OUTPUT_PATH $PROPORTIONS_PATH
echo "Annotation finished!"

sacct -j $jid1 --format=JobID,JobName,MaxRSS,Elapsed,State > $OUTPUT_PATH/preprocessing_sacct.log
sacct -j $jid2 --format=JobID,JobName,MaxRSS,Elapsed,State > $OUTPUT_PATH/sacct.log

rm -r $OUTPUT_PATH/tmp

echo "[10] FAST finished"
echo
