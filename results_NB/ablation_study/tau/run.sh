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

values=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)

mkdir -p "$OUTPUT_PATH/logs/"

for val in "${values[@]}"
do
  (
    cd /scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SNMF
    mkdir -p "$OUTPUT_PATH/v$val/"
    bash run.sh \
        "$DATA_PATH" \
        "$OUTPUT_PATH/v$val/" \
        "$val" \
        "KL_NB" \
        "$K" \
        "$PROPORTIONS_PATH" \
        "42"
  ) > "$OUTPUT_PATH/logs/SNMF_v${val}.log" 2>&1 &

  sleep 10
done

wait

source /scratch/lalonsoeste/PhD/NMF_deconvolution/.venv/bin/activate
python /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/ablation_study/tau/plot_metrics.py \
  "$OUTPUT_PATH" \
  "$PROPORTIONS_PATH"

