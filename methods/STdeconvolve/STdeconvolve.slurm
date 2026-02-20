#!/bin/bash
#SBATCH --job-name=STdeconvolve
#SBATCH --output=./logs/STdeconvolve_%j.out
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00

OUTPUT_PATH=$1
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

MARKER_REF_PATH=$2
if [ -z "$MARKER_REF_PATH" ]; then
    exit 1
fi

K=$3
if [ -z "$K" ]; then
    exit 1
fi

HUNGARIAN=$4
if [ -z "$HUNGARIAN" ]; then
    exit 1
fi

SEED=$5
if [ -z "$SEED" ]; then
    exit 1
fi

# Load R
module purge
module load R/4.4.1-gfbf-2023a
module load GSL/2.7-GCCcore-12.3.0

# Run the R script
Rscript ./STdeconvolve.R $OUTPUT_PATH $MARKER_REF_PATH $K $HUNGARIAN $SEED

