#!/bin/bash
#SBATCH --job-name=subsample_data
#SBATCH --output=/scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/sequencing_depth/logs/subsample_data.out
#SBATCH --qos=test
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:10:00

DATA_PATH=$1
if [ -z "$DATA_PATH" ]; then
    exit 1
fi

val=$2
if [ -z "$val" ]; then
    exit 1
fi

OUTPUT_PATH=$3
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

source /scratch/lalonsoeste/PhD/NMF_deconvolution/.venv/bin/activate
python /scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/sequencing_depth/subsample.py \
    $DATA_PATH \
    $val \
    "$OUTPUT_PATH/tmp/v$val.csv"