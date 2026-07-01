#!/bin/bash
#SBATCH --job-name=SNMF_cpu
#SBATCH --output=/scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/cpu/logs/SNMF_%j.out
#SBATCH --qos=regular
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=01:00:00

OUTPUT_PATH=$1
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

LOSS_FUNC=$2
if [ -z "$LOSS_FUNC" ]; then
    LOSS_FUNC="KL_NB"
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
module load CUDA
module load Apptainer

# Run the R script
apptainer exec \
    --nv \
    --bind $OUTPUT_PATH:$OUTPUT_PATH \
    --bind /scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SNMF/:/scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SNMF/ \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SNMF/r-ver_4.5.1.sif \
    Rscript /scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SNMF/SNMF.R \
    $OUTPUT_PATH \
    $LOSS_FUNC \
    $K \
    $SEED
