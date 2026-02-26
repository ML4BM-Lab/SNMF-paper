#!/bin/bash
#SBATCH --job-name=spicemix
#SBATCH --output=./logs/spicemix_%j.out
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --time=07:00:00

OUTPUT_PATH=$1
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

NITER=$2
if [ -z "$NITER" ]; then
    exit 1
fi

SEED=$3
if [ -z "$SEED" ]; then
    exit 1
fi

# Load environment
source /scratch/lalonsoeste/PhD/NMF_deconvolution/.venv/bin/activate

# Run the R script
python ./run_spicemix.py $OUTPUT_PATH $NITER $SEED

