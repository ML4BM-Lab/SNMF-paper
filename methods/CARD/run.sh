#!/bin/bash

echo "[1] CARD:"

DATA_PATH=$1
if [ -z "$DATA_PATH" ]; then
    exit 1
fi

MARKERS_PATH=$2
if [ -z "$MARKERS_PATH" ]; then
    exit 1
fi

OUTPUT_PATH=$3
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

SEED=$4
if [ -z "$SEED" ]; then
    exit 1
fi

# Load R
module purge
module load R/4.4.1-gfbf-2023a

# CARD dependencies
module load UDUNITS/2.2.28-GCCcore-12.3.0
module load GDAL/3.7.1-foss-2023a

# Make temporary folder
mkdir $OUTPUT_PATH/tmp

echo "Loading data..."
jid1=$(sbatch --parsable --wait ./load_data.slurm $DATA_PATH $MARKERS_PATH $OUTPUT_PATH)
echo "Data loaded!"

echo "Deconvolution started..."
jid2=$(sbatch --parsable --wait ./CARD.slurm $OUTPUT_PATH $SEED)
echo "Deconvolution finished!"

echo "Annotating cell types..."
Rscript ./annotation.R $DATA_PATH $MARKERS_PATH $OUTPUT_PATH
echo "Annotation finished!"

sacct -j $jid1 --format=JobID,JobName,MaxRSS,Elapsed,State > $OUTPUT_PATH/preprocessing_sacct.log
sacct -j $jid2 --format=JobID,JobName,MaxRSS,Elapsed,State > $OUTPUT_PATH/sacct.log

rm -r $OUTPUT_PATH/tmp

echo "[1] CARD finished"
echo