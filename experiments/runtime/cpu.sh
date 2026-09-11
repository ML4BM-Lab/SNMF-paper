#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

echo "[7] Spatial NMF:"

DATA_PATH="${1:-}"
if [ -z "$DATA_PATH" ]; then
    exit 1
fi

OUTPUT_PATH="${2:-}"
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

TAU="${3:-}"
if [ -z "$TAU" ]; then
    exit 1
fi

LOSS_FUNC="${4:-}"
if [ -z "$LOSS_FUNC" ]; then
    LOSS_FUNC="NB"
fi

K="${5:-}"
if [ -z "$K" ]; then
    exit 1
fi

PROPORTIONS_PATH="${6:-}"

SEED="${7:-}"
if [ -z "$SEED" ]; then
    exit 1
fi

# Load R
module purge
module load R/4.4.1-gfbf-2023a

# Make temporary folder
mkdir $OUTPUT_PATH/tmp
mkdir -p "$SCRIPT_DIR/outputs/logs"

MAX_TEST_JOBS=2
USER_NAME=$(whoami)

while true; do
    TEST_JOBS=$(squeue -u "$USER_NAME" -h -q test | wc -l)
    if (( TEST_JOBS < MAX_TEST_JOBS )); then
      echo "Loading data..."
      cd "$REPO_ROOT/methods/SNMF"
      jid1=$(sbatch --parsable --wait ./load_data.slurm  $DATA_PATH $OUTPUT_PATH $TAU)
      echo "Data loaded!"
      break
    else
      echo "Max SLURM test QoS jobs reached. Will try again in 60 seconds"
      sleep 60
    fi
done

while true; do
    TEST_JOBS=$(squeue -u "$USER_NAME" -h -q test | wc -l)
    if (( TEST_JOBS < MAX_TEST_JOBS )); then
      echo "Deconvolution started..."
      cd "$SCRIPT_DIR"
      jid2=$(sbatch --parsable --wait ./SNMF.slurm $OUTPUT_PATH $LOSS_FUNC $K $SEED)
      echo "Deconvolution finished!"
      break
    else
      echo "Max SLURM test QoS jobs reached. Will try again in 60 seconds"
      sleep 60
    fi
done

sacct -j $jid1 --format=JobID,JobName,MaxRSS,Elapsed,State > $OUTPUT_PATH/preprocessing_sacct.log
sacct -j $jid2 --format=JobID,JobName,MaxRSS,Elapsed,State > $OUTPUT_PATH/sacct.log

rm -r $OUTPUT_PATH/tmp

echo "[7] Spatial NMF finished"
echo
