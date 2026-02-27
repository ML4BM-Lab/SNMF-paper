#!/bin/bash

DATA_PATH=$1
if [ -z "$DATA_PATH" ]; then
    exit 1
fi

OUTPUT_BASE=$2
if [ -z "$OUTPUT_BASE" ]; then
    exit 1
fi

K=$3
if [ -z "$K" ]; then
    exit 1
fi

PROPORTIONS_PATH=$4
if [ -z "$PROPORTIONS_PATH" ]; then
    exit 1
fi

NRUN_VALUES=(1 2 3 4 5 6 7 8 9 10)

module load R/4.4.1-gfbf-2023a

for NRUN in "${NRUN_VALUES[@]}"; do

    echo "=== Starting run with NRUN=$NRUN ==="

    cd /scratch/lalonsoeste/PhD/NMF_deconvolution

    OUTPUT_PATH="${OUTPUT_BASE}/run_${NRUN}"
    mkdir -p "$OUTPUT_PATH"

    mkdir -p "$OUTPUT_PATH/NMF/"
    echo "NMF started..."

    jid1=$(sbatch --parsable --wait /scratch/lalonsoeste/PhD/NMF_deconvolution/experiments/starfysh_synthetic/spatial_filt_contribution/nmf.slurm $DATA_PATH $OUTPUT_PATH/NMF/ $K $NRUN)

    echo "Hungarian algorithm started..."
    Rscript /scratch/lalonsoeste/PhD/NMF_deconvolution/annotation/hungarian.R $OUTPUT_PATH/NMF/ $OUTPUT_PATH/NMF/proportions.csv $PROPORTIONS_PATH
    echo "Hungarian algorithm finished!"

    sacct -j $jid1 --format=JobID,JobName,MaxRSS,Elapsed,State > $OUTPUT_PATH/NMF/sacct.log
    
    echo "NMF finished!"

    mkdir -p "$OUTPUT_PATH/SNMF/"
    cd /scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SNMF
    bash run.sh \
        "$DATA_PATH" \
        "$OUTPUT_PATH/SNMF/" \
        $K \
        $NRUN \
        0.75 \
        "$PROPORTIONS_PATH" \
        42
    
    echo "=== Finished NRUN=$NRUN ==="
    echo ""

done

source /scratch/lalonsoeste/PhD/NMF_deconvolution/.venv/bin/activate
python /scratch/lalonsoeste/PhD/NMF_deconvolution/experiments/starfysh_synthetic/spatial_filt_contribution/results/plot.py $OUTPUT_BASE $PROPORTIONS_PATH