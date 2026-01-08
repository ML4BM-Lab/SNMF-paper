#!/bin/bash
#SBATCH --job-name=SNMF_v2
#SBATCH --output=./logs/SNMF_v2_%j.out
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00

OUTPUT_PATH=$1
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

K=$2
if [ -z "$K" ]; then
    exit 1
fi

THETA=$3
if [ -z "$THETA" ]; then
    exit 1
fi

PROBS=$4
if [ -z "$PROBS" ]; then
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
apptainer exec --nv --bind $OUTPUT_PATH:$OUTPUT_PATH ./r-ver_4.5.1.sif Rscript ./SNMF.R $OUTPUT_PATH $K $THETA $PROBS $SEED
