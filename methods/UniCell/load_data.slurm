#!/bin/bash
#SBATCH --job-name=load_data_ucd
#SBATCH --output=./logs/load_data_ucd_%j.out
#SBATCH --qos=test
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:10:00

DATA_PATH=$1
if [ -z "$DATA_PATH" ]; then
    exit 1
fi

OUTPUT_PATH=$2
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

# Load environment
source /scratch/lalonsoeste/PhD/NMF_deconvolution/.venv/bin/activate

# Run the R script
python ./load_data.py $DATA_PATH $OUTPUT_PATH