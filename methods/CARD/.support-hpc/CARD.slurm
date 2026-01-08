#!/bin/bash
#SBATCH --job-name=CARD
#SBATCH --output=./logs/CARD_%j.out
#SBATCH --qos=test
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:05:00

OUTPUT_PATH=$1
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

# Load R
module purge
module load R/4.4.1-gfbf-2023a

# CARD dependencies
module load UDUNITS/2.2.28-GCCcore-12.3.0
module load GDAL/3.7.1-foss-2023a

# Run the R script
Rscript ./CARD.R $OUTPUT_PATH

