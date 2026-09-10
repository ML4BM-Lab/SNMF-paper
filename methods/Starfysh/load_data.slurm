#!/bin/bash
#SBATCH --job-name=load_data_starfysh
#SBATCH --output=./logs/load_data_starfysh_%j.out
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:20:00

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

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

# Load environment
if [ -d "$REPO_ROOT/.venv" ]; then
    source "$REPO_ROOT/.venv/bin/activate"
fi

# Run the R script
python ./load_data.py $DATA_PATH $MARKERS_PATH $OUTPUT_PATH
