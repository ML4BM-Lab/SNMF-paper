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

echo "[9] SpiceMix"

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

NITER="${4:-}"
if [ -z "$NITER" ]; then
    exit 1
fi

SEED="${5:-}"
if [ -z "$SEED" ]; then
    exit 1
fi

DATA_PATH="$(repo_abs_path "$DATA_PATH")"
OUTPUT_PATH="$(repo_abs_path "$OUTPUT_PATH")"

cd "$SCRIPT_DIR"
mkdir -p "$SCRIPT_DIR/logs"

# Make temporary folder
mkdir -p "$OUTPUT_PATH/tmp"

MAX_TEST_JOBS=2
USER_NAME=$(whoami)

while true; do
    TEST_JOBS=$(squeue -u "$USER_NAME" -h -q test | wc -l)
    if (( TEST_JOBS < MAX_TEST_JOBS )); then
      echo "Loading data..."
      jid1=$(sbatch --parsable --wait ./load_data.slurm $DATA_PATH $OUTPUT_PATH $K)
      echo "Data loaded!"
      break
    else
      echo "Max SLURM test QoS jobs reached. Will try again in 30 seconds"
      sleep 30
    fi
done

echo "Deconvolution started..."
jid2=$(sbatch --parsable --wait ./spicemix.slurm $OUTPUT_PATH $NITER $SEED)
echo "Deconvolution finished!"

sacct -j $jid1 --format=JobID,JobName,MaxRSS,Elapsed,State > $OUTPUT_PATH/preprocessing_sacct.log
sacct -j $jid2 --format=JobID,JobName,MaxRSS,Elapsed,State > $OUTPUT_PATH/sacct.log

rm -r $OUTPUT_PATH/tmp

echo "[9] SpiceMix finished"
echo
