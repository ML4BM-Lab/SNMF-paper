#!/bin/bash
#SBATCH --job-name=mustang_step1
#SBATCH --output=mustang_step1_%j.out
#SBATCH --error=mustang_step1_%j.err
#SBATCH --partition=standard
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00

DATA_PATH=$1
if [ -z "$DATA_PATH" ]; then
    exit 1
fi

OUTPUT_PATH=$2
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

# Load R
module purge
module load R/4.4.1-gfbf-2023a

# Run the R script
Rscript ./load_data.R $DATA_PATH $OUTPUT_PATH
