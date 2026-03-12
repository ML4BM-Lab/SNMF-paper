#!/bin/bash
#SBATCH --job-name=SNMF
#SBATCH --output=./logs/SNMF_%j.out
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --constraint=rtx3090
#SBATCH --time=00:30:00

OUTPUT_PATH=$1
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

LOSS_FUNC=$2
if [ -z "$LOSS_FUNC" ]; then
    LOSS_FUNC="KL_poisson"
fi

K=$3
if [ -z "$K" ]; then
    exit 1
fi

SEED=$4
if [ -z "$SEED" ]; then
    exit 1
fi

# Load R
module purge
module load CUDA
module load Apptainer

# Run the R script
apptainer exec --nv --bind $OUTPUT_PATH:$OUTPUT_PATH r-ver_4.5.1.sif Rscript ./SNMF.R $OUTPUT_PATH $LOSS_FUNC $K $SEED
