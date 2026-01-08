#!/bin/bash
#SBATCH --job-name=ucd
#SBATCH --output=./logs/ucd%j.out
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00

OUTPUT_PATH=$1
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

SEED=$2
if [ -z "$SEED" ]; then
    exit 1
fi

# Load environment
source /scratch/lalonsoeste/PhD/NMF_deconvolution/.venv/bin/activate

# Run the R script
python ./run_ucdeconvolve.py $OUTPUT_PATH $SEED

