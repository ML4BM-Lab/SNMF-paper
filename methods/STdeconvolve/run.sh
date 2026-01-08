#!/bin/bash

echo "[5] STdeconvolve:"

DATA_PATH=$1
if [ -z "$DATA_PATH" ]; then
    exit 1
fi

MARKER_REF_PATH=$2
if [ -z "$MARKER_REF_PATH" ]; then
    exit 1
fi

OUTPUT_PATH=$3
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

K=$4
if [ -z "$K" ]; then
    exit 1
fi

HUNGARIAN=$5
if [ -z "$HUNGARIAN" ]; then
    exit 1
fi

SEED=$6
if [ -z "$SEED" ]; then
    exit 1
fi

# Load R
module purge
module load R/4.4.1-gfbf-2023a
module load GSL/2.7-GCCcore-12.3.0

# Make temporary folder
mkdir $OUTPUT_PATH/tmp

MAX_TEST_JOBS=2
USER_NAME=$(whoami)

while true; do
    TEST_JOBS=$(squeue -u "$USER_NAME" -h -q test | wc -l)
    if (( TEST_JOBS < MAX_TEST_JOBS )); then
      echo "Loading data..."
      jid1=$(sbatch --parsable --wait ./load_data.slurm $DATA_PATH $OUTPUT_PATH)
      echo "Data loaded!"
      break
    else
      echo "Max SLURM test QoS jobs reached. Will try again in 30 seconds"
      sleep 30
    fi
done

echo "Deconvolution started..."
jid2=$(sbatch --parsable --wait ./STdeconvolve.slurm $OUTPUT_PATH $MARKER_REF_PATH $K $HUNGARIAN $SEED)
echo "Deconvolution finished!"

sacct -j $jid1 --format=JobID,JobName,MaxRSS,Elapsed,State > $OUTPUT_PATH/preprocessing_sacct.log
sacct -j $jid2 --format=JobID,JobName,MaxRSS,Elapsed,State > $OUTPUT_PATH/sacct.log

rm -r $OUTPUT_PATH/tmp

echo "[5] STdeconvolve finished"
echo