#!/bin/bash
#SBATCH --job-name=bayestme
#SBATCH --output=./logs/bayestme_%j.out
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

OUTPUT_PATH=$1
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

K=$2
if [ -z "$K" ]; then
    exit 1
fi

RHO=$3
if [ -z "$RHO" ]; then
    exit 1
fi

SEED=$4
if [ -z "$SEED" ]; then
    exit 1
fi

# Load environment
if [ -d "$REPO_ROOT/.venv" ]; then
    source "$REPO_ROOT/.venv/bin/activate"
fi

# Run the R script
python ./run_bayestme.py $OUTPUT_PATH $K $RHO $SEED
