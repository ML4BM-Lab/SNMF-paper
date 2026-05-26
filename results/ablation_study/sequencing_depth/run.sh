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

mkdir -p "$OUTPUT_PATH/tmp"
mkdir -p "$OUTPUT_PATH/logs"

for val in "${values[@]}"
do
  mkdir -p "$OUTPUT_PATH/v$val"

  echo "Subsampling data..."
  sbatch --parsable --wait /scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/sequencing_depth/subsample_data.slurm \
    $DATA_PATH \
    $val \
    "$OUTPUT_PATH"
  echo "Data subsampled!"

  (
    cd /scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SNMF
    bash run.sh \
        "$OUTPUT_PATH/tmp/v$val.csv" \
        "$OUTPUT_PATH/v$val/" \
        "0.4" \
        "KL_poisson" \
        "$K" \
        "$PROPORTIONS_PATH" \
        "42"
  ) > "$OUTPUT_PATH/logs/SNMF_v${val}.log" 2>&1 &
  echo "SNMF Launched"

  sleep 10
done

wait

rm -r "$OUTPUT_PATH/tmp"

source /scratch/lalonsoeste/PhD/NMF_deconvolution/.venv/bin/activate
python /scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/tau/plot_metrics.py \
  "$OUTPUT_PATH" \
  "$PROPORTIONS_PATH"

