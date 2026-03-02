#!/usr/bin/env bash
set -euo pipefail

# Default values
DATA_PATH=""
MARKERS_PATH=""
OUTPUT_PATH=""
VISIUM=false
K=""
PROPORTIONS_PATH=""
SNMF_VALUE="0.5"
STARFYSH_LR="1e-6"
HUNGARIAN=false
SEED=42

# Parse flags
while [[ $# -gt 0 ]]; do
  case "$1" in
    --data_path=*)
      DATA_PATH="${1#*=}"
      ;;
    --markers_path=*)
      MARKERS_PATH="${1#*=}"
      ;;
    --output_path=*)
      OUTPUT_PATH="${1#*=}"
      ;;
    --visium=*)
      VISIUM="${1#*=}"
      ;;
    --k=*)
      K="${1#*=}"
      ;;
    --proportions_path=*)
      PROPORTIONS_PATH="${1#*=}"
      ;;
    --snmf_value=*)
      SNMF_VALUE="${1#*=}"
      ;;
    --hungarian=*)
      HUNGARIAN="${1#*=}"
      ;;
    --starfysh_lr=*)
      STARFYSH_LR="${1#*=}"
      ;;
    --seed=*)
      SEED="${1#*=}"
      ;;
    --help|-h)
      echo "Usage: $0 --data_path=FILE --markers_path=FILE --output_path=DIR --k=INT [--proportions_path=FILE] [--visium=true|false] --starfysh_lr=STARFYSH_LR] [--hungarian=true|false] [-seed=SEED]"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
  shift
done

# Validate required args
if [ -z "$DATA_PATH" ] || [ -z "$MARKERS_PATH" ] || [ -z "$OUTPUT_PATH" ] || [ -z "$K" ]; then
  echo "Error: Missing required arguments."
  echo "Usage: $0 --data_path=FILE --markers_path=FILE --output_path=DIR --k=INT [--starfysh_lr=STARFYSH_LR] [--hungarian=true|false] [-ssed=SEED]"
  exit 1
fi

# Pipeline
cd "/scratch/lalonsoeste/PhD/NMF_deconvolution/"
mkdir -p "$OUTPUT_PATH/logs"

## RETROFIT
(
  cd /scratch/lalonsoeste/PhD/NMF_deconvolution/methods/RETROFIT
  mkdir -p "$OUTPUT_PATH/RETROFIT/"
  bash run.sh \
      "$DATA_PATH" \
      "$K" \
      "$MARKERS_PATH" \
      "$OUTPUT_PATH/RETROFIT/" \
      $SEED
) > "$OUTPUT_PATH/logs/retrofit.log" 2>&1 &

sleep 10 # This avoids 'sbatch: error: Batch job submission failed: Socket timed out on send/recv operation'

## SNMF (ours)
(
  cd /scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SNMF
  mkdir -p "$OUTPUT_PATH/SNMF/"
  bash run.sh \
      "$DATA_PATH" \
      "$OUTPUT_PATH/SNMF/" \
      "$SNMF_VALUE" \
      $K \
      10 \
      0.75 \
      "$PROPORTIONS_PATH" \
      $SEED
) > "$OUTPUT_PATH/logs/SNMF.log" 2>&1 &

sleep 10

## STdeconvolve 
(
  cd /scratch/lalonsoeste/PhD/NMF_deconvolution/methods/STdeconvolve
  mkdir -p "$OUTPUT_PATH/STdeconvolve/"
  bash run.sh \
      "$DATA_PATH" \
      "$MARKERS_PATH" \
      "$OUTPUT_PATH/STdeconvolve/" \
      $K \
      true \
      $SEED
) > "$OUTPUT_PATH/logs/STdeconvolve.log" 2>&1 &

sleep 10

## SMART
(
  cd /scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SMART
  mkdir -p "$OUTPUT_PATH/SMART/"
  bash run.sh \
      "$DATA_PATH" \
      "$MARKERS_PATH" \
      "$OUTPUT_PATH/SMART/" \
      $SEED
) > "$OUTPUT_PATH/logs/SMART.log" 2>&1 &

sleep 10

## starfysh
(
  cd /scratch/lalonsoeste/PhD/NMF_deconvolution/methods/Starfysh
  mkdir -p "$OUTPUT_PATH/starfysh/"
  bash run.sh \
      "$DATA_PATH" \
      "$MARKERS_PATH" \
      "$OUTPUT_PATH/starfysh/" \
      $STARFYSH_LR \
      $SEED
) > "$OUTPUT_PATH/logs/starfysh.log" 2>&1 &

sleep 10

## BayesTME
(
  cd /scratch/lalonsoeste/PhD/NMF_deconvolution/methods/BayesTME
  mkdir -p "$OUTPUT_PATH/BayesTME/"
  bash run.sh \
      "$DATA_PATH" \
      "$OUTPUT_PATH/BayesTME/" \
      "$VISIUM" \
      $K \
      0.5 \
      $SEED
) > "$OUTPUT_PATH/logs/bayestme.log" 2>&1 &

sleep 10

## SpiceMix
(
  cd /scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SpiceMix
  mkdir -p "$OUTPUT_PATH/SpiceMix/"
  bash run.sh \
      "$DATA_PATH" \
      "$OUTPUT_PATH/SpiceMix/" \
      $K \
      200 \
      $SEED
) > "$OUTPUT_PATH/logs/spicemix.log" 2>&1 &

sleep 10

## CARD
(
  cd /scratch/lalonsoeste/PhD/NMF_deconvolution/methods/CARD
  mkdir -p "$OUTPUT_PATH/CARD/"
  bash run.sh \
      "$DATA_PATH" \
      "$MARKERS_PATH" \
      "$OUTPUT_PATH/CARD/" \
      $SEED
) > "$OUTPUT_PATH/logs/CARD.log" 2>&1 &


# wait for all background jobs to finish
wait

# Hungarian annotation
if [[ "$HUNGARIAN" == "true" ]]; then
  module load R/4.4.1-gfbf-2023a
  for method in CARD RETROFIT STdeconvolve SMART starfysh BayesTME SpiceMix SNMF; do
    if [[ -f "$OUTPUT_PATH/$method/${method}_proportions.csv" ]]; then
      echo "Computing hungarian algorithm for $method ..."
      Rscript ./annotation/hungarian.R \
          "$OUTPUT_PATH/$method/" \
          "$OUTPUT_PATH/$method/${method}_proportions.csv" \
          "$PROPORTIONS_PATH"
    fi
  done
fi

# Plot results
if [ ! -z "$PROPORTIONS_PATH" ]; then
  source /scratch/lalonsoeste/PhD/NMF_deconvolution/.venv/bin/activate
  mkdir -p "$OUTPUT_PATH/plots"
  python /scratch/lalonsoeste/PhD/NMF_deconvolution/analysis/plot_metrics.py \
      "$OUTPUT_PATH" \
      "$PROPORTIONS_PATH" \
      "$HUNGARIAN"

  python /scratch/lalonsoeste/PhD/NMF_deconvolution/analysis/plot_proportions.py \
      "$OUTPUT_PATH" \
      "$PROPORTIONS_PATH" \
      "$HUNGARIAN"
fi
