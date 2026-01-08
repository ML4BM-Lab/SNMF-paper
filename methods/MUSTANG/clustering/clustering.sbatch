#!/bin/bash
#SBATCH --job-name=mustang_step3
#SBATCH --output=mustang_step3_%j.out
#SBATCH --error=mustang_step3_%j.err
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00

OUTPUT_PATH=$1
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

module purge
module load R/4.4.1-gfbf-2023a

Rscript ./clustering.R $OUTPUT_PATH
