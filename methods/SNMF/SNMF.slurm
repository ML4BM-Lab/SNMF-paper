#!/bin/bash
#SBATCH --job-name=SNMF
#SBATCH --output=./logs/SNMF_%j.out
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --gres=gpu:1
#SBATCH --constraint=rtx3090
#SBATCH --time=03:00:00

OUTPUT_PATH=$1
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

LOSS_FUNC=$2
if [ -z "$LOSS_FUNC" ]; then
    LOSS_FUNC="KL_poisson"
fi

DISPERSION_MODE=$3
if [ -z "$DISPERSION_MODE" ]; then
    DISPERSION_MODE="full"
fi

K=$4
if [ -z "$K" ]; then
    exit 1
fi

SEED=$5
if [ -z "$SEED" ]; then
    exit 1
fi

# Load R
module purge
module load CUDA
module load Apptainer

# Run the R script
apptainer exec --nv --bind $OUTPUT_PATH:$OUTPUT_PATH r-ver_4.5.1.sif Rscript ./SNMF.R $OUTPUT_PATH $LOSS_FUNC $DISPERSION_MODE $K $SEED
