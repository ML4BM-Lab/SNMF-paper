#!/bin/bash
#SBATCH --job-name=ucd
#SBATCH --output=./logs/ucd%j.out
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --constraint=rtx3090
#SBATCH --time=01:00:00

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

OUTPUT_PATH=$1
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

SEED=$2
if [ -z "$SEED" ]; then
    exit 1
fi

# Load environment
if [ -d "$REPO_ROOT/.venv" ]; then
    source "$REPO_ROOT/.venv/bin/activate"
fi

# Run the R script
python ./run_ucdeconvolve.py $OUTPUT_PATH $SEED
