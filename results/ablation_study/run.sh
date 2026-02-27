#/bin/bash

DATA_PATH="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/counts.csv"
OUTPUT_PATH="/scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/TNBC"
K="5"
PROPORTIONS_PATH="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/proportions.csv"
SEED=42

gammas=(1 2 3 4 5 6 7)

for gamma in "${gammas[@]}"
do
  (
    cd /scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SNMF
    mkdir -p "$OUTPUT_PATH/gamma_$gamma/"
    bash run.sh \
        "$DATA_PATH" \
        "$OUTPUT_PATH/gamma_$gamma/" \
        "$gamma" \
        "$K" \
        10 \
        0.75 \
        "$PROPORTIONS_PATH" \
        "$SEED"
  ) > "$OUTPUT_PATH/logs/SNMF_gamma_${gamma}.log" 2>&1 &

  sleep 10
done

