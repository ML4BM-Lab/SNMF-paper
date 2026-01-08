#!/bin/bash
#SBATCH --job-name=load_data_STdeconvolve
#SBATCH --output=./logs/load_data_STdeconvolve_%j.out
#SBATCH --qos=test
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:10:00

DATA_PATH=$1
if [ -z "$DATA_PATH" ]; then
    exit 1
fi

OUTPUT_PATH=$2
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

# Load environment
module purge
module load R/4.4.1-gfbf-2023a
module load GSL/2.7-GCCcore-12.3.0

# Run the R script
Rscript ./load_data.R $DATA_PATH $OUTPUT_PATH