#!/bin/bash
#SBATCH --job-name=starfysh
#SBATCH --output=./logs/starfysh_%j.out
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00

OUTPUT_PATH=$1
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

LEARNING_RATE=$2
if [ -z "$LEARNING_RATE" ]; then
    exit 1
fi

SEED=$3
if [ -z "$SEED" ]; then
    exit 1
fi

# Load environment
source /scratch/lalonsoeste/PhD/NMF_deconvolution/.venv/bin/activate

# Run the python script
python ./run_starfysh.py $OUTPUT_PATH $LEARNING_RATE $SEED

