#!/bin/bash
#SBATCH --job-name=subsample_data
#SBATCH --output=/scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/ablation_study/ngenes/logs/subsample_data.out
#SBATCH --qos=test
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:10:00

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

source /scratch/lalonsoeste/PhD/NMF_deconvolution/.venv/bin/activate
python /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/ablation_study/ngenes/subsample.py \
    $DATA_PATH \
    $MARKERS_PATH \
    $val \
    "$OUTPUT_PATH/tmp/v$val.csv" \
    "$OUTPUT_PATH/tmp/mg$val.csv"