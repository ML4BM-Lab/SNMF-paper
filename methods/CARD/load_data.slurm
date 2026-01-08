#!/bin/bash
#SBATCH --job-name=load_data_CARD
#SBATCH --output=./logs/load_data_CARD_%j.out
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --time=00:05:00

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

# Load R
module purge
module load R/4.4.1-gfbf-2023a

# CARD dependencies
module load UDUNITS/2.2.28-GCCcore-12.3.0
module load GDAL/3.7.1-foss-2023a

# Run the R script
Rscript ./load_data.R $DATA_PATH $MARKERS_PATH $OUTPUT_PATH