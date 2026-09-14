#!/bin/bash
#SBATCH --job-name=FAST
#SBATCH --output=./logs/FAST_%j.out
#SBATCH --qos=test
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:10:00

OUTPUT_PATH=$1
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

SEED=$2
if [ -z "$SEED" ]; then
    exit 1
fi

# Load R
module purge
module load R/4.4.1-gfbf-2023a

# Run the R script
Rscript ./FAST.R $OUTPUT_PATH $SEED

