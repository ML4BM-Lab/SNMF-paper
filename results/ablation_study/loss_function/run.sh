#/bin/bash

# DATA_PATH="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/counts.csv"
# OUTPUT_PATH="/scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/loss_function/TNBC"
# K="5"
# PROPORTIONS_PATH="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/proportions.csv"

DATA_PATH="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/mixture_file.csv"
OUTPUT_PATH="/scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/loss_function/PDAC"
K="20"
PROPORTIONS_PATH="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/proportions.csv"

SEED=42

loss_functions=(KL_poisson squared_error KL_NB)

for loss_func in "${loss_functions[@]}"
do
  (
    cd /scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SNMF
    mkdir -p "$OUTPUT_PATH/$loss_func/"
    bash run.sh \
        "$DATA_PATH" \
        "$OUTPUT_PATH/$loss_func/" \
        "0.5" \
        "$loss_func" \
        "$K" \
        "$PROPORTIONS_PATH" \
        "$SEED"
  ) > "$OUTPUT_PATH/logs/SNMF_${loss_func}.log" 2>&1 &
done

wait

# source /scratch/lalonsoeste/PhD/NMF_deconvolution/.venv/bin/activate
# python /scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/tau/plot_metrics.py \
#   "$OUTPUT_PATH" \
#   "$PROPORTIONS_PATH"

