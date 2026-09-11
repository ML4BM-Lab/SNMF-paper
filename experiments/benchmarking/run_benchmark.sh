#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
VENV_DIR="${REPO_ROOT}/.venv"
SCRIPTS_DIR="${SCRIPT_DIR}/scripts"
METHODS_DIR="${REPO_ROOT}/methods"

repo_abs_path() {
  case "$1" in
    /*) printf "%s\n" "$1" ;;
    *) printf "%s/%s\n" "$REPO_ROOT" "$1" ;;
  esac
}

# Default values
DATA_PATH=""
MARKERS_PATH=""
OUTPUT_PATH=""
VISIUM=false
K=""
PROPORTIONS_PATH=""
SNMF_TAU="0.8"
STARFYSH_LR="1e-6"
HUNGARIAN=true
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
    --snmf_tau=*)
      SNMF_TAU="${1#*=}"
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
      echo "Usage: $0 --data_path=FILE --markers_path=FILE --output_path=DIR --k=INT [--snmf_tau=SNMF_TAU] [--proportions_path=FILE] [--visium=true|false] [--starfysh_lr=STARFYSH_LR] [--hungarian=true|false] [--seed=SEED]"
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
  echo "Usage: $0 --data_path=FILE --markers_path=FILE --output_path=DIR --k=INT [--proportions_path=FILE] [--visium=true|false] [--starfysh_lr=STARFYSH_LR] [--hungarian=true|false] [--seed=SEED]"
  exit 1
fi

DATA_PATH="$(repo_abs_path "$DATA_PATH")"
MARKERS_PATH="$(repo_abs_path "$MARKERS_PATH")"
OUTPUT_PATH="$(repo_abs_path "$OUTPUT_PATH")"
if [ -n "$PROPORTIONS_PATH" ]; then
  PROPORTIONS_PATH="$(repo_abs_path "$PROPORTIONS_PATH")"
fi

# Pipeline
cd "$REPO_ROOT"
mkdir -p "$OUTPUT_PATH/logs"

## RETROFIT
(
  cd "$METHODS_DIR/RETROFIT"
  mkdir -p "$OUTPUT_PATH/RETROFIT/"
  bash run.sh \
      "$DATA_PATH" \
      "$K" \
      "$MARKERS_PATH" \
      "$OUTPUT_PATH/RETROFIT/" \
      $SEED
) > "$OUTPUT_PATH/logs/retrofit.log" 2>&1 &

sleep 10

## SNMF (ours)
(
  cd "$METHODS_DIR/SNMF"
  mkdir -p "$OUTPUT_PATH/SNMF/"
  bash run.sh \
      "$DATA_PATH" \
      "$OUTPUT_PATH/SNMF/" \
      "$SNMF_TAU" \
      "NB" \
      $K \
      "$PROPORTIONS_PATH" \
      $SEED \
      "full"
) > "$OUTPUT_PATH/logs/SNMF.log" 2>&1 &

sleep 10

## NMF (S=I)
(
  cd "$METHODS_DIR/SNMF"
  mkdir -p "$OUTPUT_PATH/NMF/"
  bash run.sh \
      "$DATA_PATH" \
      "$OUTPUT_PATH/NMF/" \
      1 \
      "NB" \
      $K \
      "$PROPORTIONS_PATH" \
      $SEED \
      "full"
) > "$OUTPUT_PATH/logs/NMF.log" 2>&1 &

sleep 10

## STdeconvolve 
(
  cd "$METHODS_DIR/STdeconvolve"
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
  cd "$METHODS_DIR/SMART"
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
  cd "$METHODS_DIR/Starfysh"
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
  cd "$METHODS_DIR/BayesTME"
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
  cd "$METHODS_DIR/SpiceMix"
  mkdir -p "$OUTPUT_PATH/SpiceMix/"
  bash run.sh \
      "$DATA_PATH" \
      "$OUTPUT_PATH/SpiceMix/" \
      $K \
      100 \
      $SEED
) > "$OUTPUT_PATH/logs/spicemix.log" 2>&1 &

sleep 10

## CARD
(
  cd "$METHODS_DIR/CARD"
  mkdir -p "$OUTPUT_PATH/CARD/"
  bash run.sh \
      "$DATA_PATH" \
      "$MARKERS_PATH" \
      "$OUTPUT_PATH/CARD/" \
      $SEED
) > "$OUTPUT_PATH/logs/CARD.log" 2>&1 &


# wait for all background jobs to finish
wait

mv "$OUTPUT_PATH/NMF/SNMF_proportions.csv" "$OUTPUT_PATH/NMF/NMF_proportions.csv"

# Hungarian annotation
if [[ "$HUNGARIAN" == "true" ]]; then
  module load R/4.4.1-gfbf-2023a
  for method in CARD RETROFIT STdeconvolve SMART starfysh BayesTME SpiceMix NMF SNMF; do
    if [[ -f "$OUTPUT_PATH/$method/${method}_proportions.csv" ]]; then
      echo "Computing hungarian algorithm for $method ..."
      Rscript "$SCRIPTS_DIR/hungarian.R" \
          "$OUTPUT_PATH/$method/" \
          "$OUTPUT_PATH/$method/${method}_proportions.csv" \
          "$PROPORTIONS_PATH"
    fi
  done
fi

# Plot results
if [ ! -z "$PROPORTIONS_PATH" ]; then
  if [ -d "$VENV_DIR" ]; then
    source "$VENV_DIR/bin/activate"
  fi
  mkdir -p "$OUTPUT_PATH/plots"
  python "$SCRIPTS_DIR/plot_metrics.py" \
      "$OUTPUT_PATH" \
      "$PROPORTIONS_PATH" \
      "$HUNGARIAN"

  python "$SCRIPTS_DIR/plot_proportions.py" \
      "$OUTPUT_PATH" \
      "$PROPORTIONS_PATH" \
      "$HUNGARIAN"
fi
