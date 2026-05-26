#!/bin/bash

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

seeds=(1 2 3 4 5)

mkdir -p "$OUTPUT_PATH/logs"

for seed in "${seeds[@]}"
do
  (
    cd /scratch/lalonsoeste/PhD/NMF_deconvolution/results/cpu
    mkdir -p "$OUTPUT_PATH/cpu/$seed/"
    bash cpu.sh \
        "$DATA_PATH" \
        "$OUTPUT_PATH/cpu/$seed/" \
        "0.4" \
        "KL_poisson" \
        "$K" \
        "$PROPORTIONS_PATH" \
        "$seed"
  ) > "$OUTPUT_PATH/logs/cpu_${seed}.log" 2>&1 &

  sleep 10

  (
    cd /scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SNMF
    mkdir -p "$OUTPUT_PATH/gpu/$seed/"
    bash run.sh \
        "$DATA_PATH" \
        "$OUTPUT_PATH/gpu/$seed/" \
        "0.4" \
        "KL_poisson" \
        "$K" \
        "$PROPORTIONS_PATH" \
        "$seed"
  ) > "$OUTPUT_PATH/logs/gpu_${seed}.log" 2>&1 &

  sleep 10
done

wait

source /scratch/lalonsoeste/PhD/NMF_deconvolution/.venv/bin/activate
python /scratch/lalonsoeste/PhD/NMF_deconvolution/results/cpu/plot_time.py \
  "$OUTPUT_PATH"