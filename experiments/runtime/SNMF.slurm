#!/bin/bash
#SBATCH --job-name=SNMF_cpu
#SBATCH --output=./outputs/logs/SNMF_%j.out
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=01:00:00

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SNMF_METHOD_DIR="${REPO_ROOT}/methods/SNMF"

OUTPUT_PATH=$1
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

LOSS_FUNC=$2
if [ -z "$LOSS_FUNC" ]; then
    LOSS_FUNC="NB"
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
module load Apptainer

# Run the R script
apptainer exec \
    --bind "$OUTPUT_PATH:$OUTPUT_PATH" \
    --bind "$SNMF_METHOD_DIR:$SNMF_METHOD_DIR" \
    "$SNMF_METHOD_DIR/r-ver_4.5.1.sif" \
    Rscript "$SNMF_METHOD_DIR/SNMF.R" \
    "$OUTPUT_PATH" \
    "$LOSS_FUNC" \
    "full" \
    "$K" \
    "$SEED"
