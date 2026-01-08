#!/bin/bash
#SBATCH --job-name=mustang_step4
#SBATCH --output=mustang_step4_%j.out
#SBATCH --error=mustang_step4_%j.err
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00

OUTPUT_PATH=$1
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

module purge
module load R/4.4.1-gfbf-2023a

Rscript transcriptional_graph.R $OUTPUT_PATH

source ../../.venv/bin/activate
python spatial_graph.py $OUTPUT_PATH

Rscript similarity_graph.R $OUTPUT_PATH
