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

loss_functions=(squared_error KL_poisson KL_NB)

mkdir -p "$OUTPUT_PATH/logs/"

for loss_func in "${loss_functions[@]}"
do
  (
    cd /scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SNMF
    mkdir -p "$OUTPUT_PATH/$loss_func/"
    bash run.sh \
        "$DATA_PATH" \
        "$OUTPUT_PATH/$loss_func/" \
        "0.4" \
        "$loss_func" \
        "$K" \
        "$PROPORTIONS_PATH" \
        "42"
  ) > "$OUTPUT_PATH/logs/SNMF_${loss_func}.log" 2>&1 &

  sleep 10
  
done

wait

source /scratch/lalonsoeste/PhD/NMF_deconvolution/.venv/bin/activate
python /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/ablation_study/loss_function/plot_metrics.py \
  "$OUTPUT_PATH" \
  "$PROPORTIONS_PATH"

