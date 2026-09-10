#!/bin/bash
#SBATCH --job-name=subsample_data
#SBATCH --output=./outputs/logs/subsample_data_%j.out
#SBATCH --qos=test
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:10:00

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
VENV_DIR="${REPO_ROOT}/.venv"

DATA_PATH=$1
if [ -z "$DATA_PATH" ]; then
    exit 1
fi

MARKERS_PATH=$2
if [ -z "$MARKERS_PATH" ]; then
    exit 1
fi

val=$3
if [ -z "$val" ]; then
    exit 1
fi

OUTPUT_PATH=$4
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

mkdir -p "$SCRIPT_DIR/outputs/logs"

if [ -d "$VENV_DIR" ]; then
    source "$VENV_DIR/bin/activate"
fi
python "$SCRIPT_DIR/subsample.py" \
    "$DATA_PATH" \
    "$MARKERS_PATH" \
    "$val" \
    "$OUTPUT_PATH/tmp/v$val.csv" \
    "$OUTPUT_PATH/tmp/mg$val.csv"
