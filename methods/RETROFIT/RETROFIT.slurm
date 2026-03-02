#!/bin/bash
#SBATCH --job-name=RETROFIT
#SBATCH --output=./logs/RETROFIT_%j.out
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00

OUTPUT_PATH=$1
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

K=$2
if [ -z "$K" ]; then
    exit 1
fi

SEED=$3
if [ -z "$SEED" ]; then
    exit 1
fi

# Load R
module purge
module load R/4.4.1-gfbf-2023a

# Run the R script
Rscript ./RETROFIT.R $OUTPUT_PATH $K $SEED

