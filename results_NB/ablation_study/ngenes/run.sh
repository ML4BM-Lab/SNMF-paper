#/bin/bash

DATA_PATH=$1
if [ -z "$DATA_PATH" ]; then
    exit 1
fi

MARKERS_PATH=$2
if [ -z "$MARKERS_PATH" ]; then
    exit 1
fi

SNMF_TAU=$3
if [ -z "$SNMF_TAU" ]; then
    exit 1
fi

OUTPUT_PATH=$4
if [ -z "$OUTPUT_PATH" ]; then
    exit 1
fi

K=$5
if [ -z "$K" ]; then
    exit 1
fi

PROPORTIONS_PATH=$6
if [ -z "$PROPORTIONS_PATH" ]; then
    exit 1
fi

ngenes=(100 200 400 800 1600 3200 6400 12800 25600)

mkdir -p "$OUTPUT_PATH/tmp"
mkdir -p "$OUTPUT_PATH/logs"

for val in "${ngenes[@]}"
do
  mkdir -p "$OUTPUT_PATH/v$val"

  echo "Subsampling genes..."
  sbatch --parsable --wait /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/ablation_study/ngenes/subsample_data.slurm \
    $DATA_PATH \
    $MARKERS_PATH \
    $val \
    "$OUTPUT_PATH"
  echo "Genes subsampled!"

  (
    bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/benchmarking/benchmark.sh \
        --data_path="$OUTPUT_PATH/tmp/v$val.csv" \
        --markers_path="$OUTPUT_PATH/tmp/mg$val.csv" \
        --output_path="$OUTPUT_PATH/v$val/" \
        --k="$K" \
        --proportions_path="$PROPORTIONS_PATH" \
        --snmf_tau="$SNMF_TAU" \
        --starfysh_lr=1e-6 \
  ) > "$OUTPUT_PATH/logs/SNMF_v${val}.log" 2>&1 &

  echo "SNMF Launched"

  sleep 100
done

wait

rm -r "$OUTPUT_PATH/tmp"

source /scratch/lalonsoeste/PhD/NMF_deconvolution/.venv/bin/activate
python /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/ablation_study/ngenes/plot_metrics.py \
  "$OUTPUT_PATH"
