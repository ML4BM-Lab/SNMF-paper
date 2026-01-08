#!/bin/bash
#SBATCH --job-name=NMF
#SBATCH --output=/scratch/lalonsoeste/PhD/NMF_deconvolution/experiments/starfysh_synthetic/spatial_filt_contribution/logs/NMF_%j.out
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=05:00:00

DATA_PATH=$1
if [ -z "$DATA_PATH" ]; then
    exit 1
fi

OUTPUT_PATH=$2
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

K=$3
if [ -z "$K" ]; then
    exit 1
fi

NRUN=$4
if [ -z "$NRUN" ]; then
    exit 1
fi

# Load R
module purge
module load R/4.4.1-gfbf-2023a

# Run the R script
Rscript /scratch/lalonsoeste/PhD/NMF_deconvolution/experiments/starfysh_synthetic/spatial_filt_contribution/nmf.R $DATA_PATH $OUTPUT_PATH $K $NRUN
