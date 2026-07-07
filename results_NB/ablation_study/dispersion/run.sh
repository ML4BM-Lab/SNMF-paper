#/bin/bash

DATA_PATH=$1
if [ -z "$DATA_PATH" ]; then
    exit 1
fi

OUTPUT_PATH=$2
if [ -z "$OUTPUT_PATH" ]; then
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

dispersion_modes=(full gene spot)

mkdir -p "$OUTPUT_PATH/logs/"

for mode in "${dispersion_modes[@]}"
do
  (
    cd /scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SNMF
    mkdir -p "$OUTPUT_PATH/$mode/"
    bash run.sh \
        "$DATA_PATH" \
        "$OUTPUT_PATH/$mode/" \
        "0.8" \
        "KL_NB" \
        "$K" \
        "$PROPORTIONS_PATH" \
        "42" \
        "$mode"
  ) > "$OUTPUT_PATH/logs/SNMF_${mode}.log" 2>&1 &

  sleep 10
  
done

wait

source /scratch/lalonsoeste/PhD/NMF_deconvolution/.venv/bin/activate
python /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/ablation_study/dispersion/plot_metrics.py \
  "$OUTPUT_PATH" \
  "$PROPORTIONS_PATH" \
  "$DATA_PATH"

