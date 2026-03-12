#/bin/bash

# DATA_PATH="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/counts.csv"
# OUTPUT_PATH="/scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/tau/TNBC"
# K="5"
# PROPORTIONS_PATH="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/proportions.csv"

DATA_PATH="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/mixture_file.csv"
OUTPUT_PATH="/scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/tau/PDAC"
K="20"
PROPORTIONS_PATH="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/proportions.csv"

SEED=42

values=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)

for val in "${values[@]}"
do
  (
    cd /scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SNMF
    mkdir -p "$OUTPUT_PATH/v$val/"
    bash run.sh \
        "$DATA_PATH" \
        "$OUTPUT_PATH/v$val/" \
        "$val" \
        "KL_poisson" \
        "$K" \
        "$PROPORTIONS_PATH" \
        "$SEED"
  ) > "$OUTPUT_PATH/logs/SNMF_v${val}.log" 2>&1 &

  sleep 10
done

wait

source /scratch/lalonsoeste/PhD/NMF_deconvolution/.venv/bin/activate
python /scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/tau/plot_metrics.py \
  "$OUTPUT_PATH" \
  "$PROPORTIONS_PATH"

