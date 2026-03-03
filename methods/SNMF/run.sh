#!/bin/bash

echo "[7] Spatial NMF:"

DATA_PATH=$1
if [ -z "$DATA_PATH" ]; then
    exit 1
fi

OUTPUT_PATH=$2
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

TAU=$3
if [ -z "$TAU" ]; then
    exit 1
fi

K=$4
if [ -z "$K" ]; then
    exit 1
fi

PROPORTIONS_PATH=$5

SEED=$6
if [ -z "$SEED" ]; then
    exit 1
fi

# Load R
module purge
module load R/4.4.1-gfbf-2023a

# Make temporary folder
mkdir $OUTPUT_PATH/tmp

MAX_TEST_JOBS=2
USER_NAME=$(whoami)

while true; do
    TEST_JOBS=$(squeue -u "$USER_NAME" -h -q test | wc -l)
    if (( TEST_JOBS < MAX_TEST_JOBS )); then
      echo "Loading data..."
      jid1=$(sbatch --parsable --wait ./load_data.slurm  $DATA_PATH $OUTPUT_PATH $TAU)
      echo "Data loaded!"
      break
    else
      echo "Max SLURM test QoS jobs reached. Will try again in 30 seconds"
      sleep 30
    fi
done

echo "Deconvolution started..."
jid2=$(sbatch --parsable --wait ./SNMF.slurm $OUTPUT_PATH $K $SEED)
echo "Deconvolution finished!"

if [ ! -z "$PROPORTIONS_PATH" ]; then
    echo "Annotation started..."
    Rscript ./annotate.R $OUTPUT_PATH $PROPORTIONS_PATH
    echo "Annotation finished!"
fi

sacct -j $jid1 --format=JobID,JobName,MaxRSS,Elapsed,State > $OUTPUT_PATH/preprocessing_sacct.log
sacct -j $jid2 --format=JobID,JobName,MaxRSS,Elapsed,State > $OUTPUT_PATH/sacct.log

rm -r $OUTPUT_PATH/tmp

echo "[7] Spatial NMF finished"
echo