#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
DATA_DIR="${REPO_ROOT}/data/DLPFC/final"
OUTPUT_DIR="${REPO_ROOT}/experiments/benchmarking/outputs/DLPFC"
VENV_DIR="${REPO_ROOT}/.venv"
PROCESS_DATA="${REPO_ROOT}/data/DLPFC/process_data.py"

if [ -d "$VENV_DIR" ]; then
    source "$VENV_DIR/bin/activate"
fi

mkdir -p "$OUTPUT_DIR"

sample_dirs=("$DATA_DIR"/*)
if [[ ${#sample_dirs[@]} -eq 0 || ! -d "${sample_dirs[0]}" ]]; then
    echo "No DLPFC sample directories found in $DATA_DIR"
    exit 1
fi

expected_k_for_sample() {
    case "$1" in
        151669|151670|151671|151672) printf "5\n" ;;
        *) printf "7\n" ;;
    esac
}

for sample_dir in "${sample_dirs[@]}"; do
    [ -d "$sample_dir" ] || continue
    sample="$(basename "$sample_dir")"
    counts_path="${sample_dir}/${sample}_counts.csv"

    [ -e "$counts_path" ] || {
        echo "No DLPFC count matrices found in $DATA_DIR"
        exit 1
    }

    markers_path="${sample_dir}/${sample}_marker_genes.csv"
    k_path="${sample_dir}/${sample}_K.txt"
    h5ad_path="${sample_dir}/${sample}_filtered_feature_bc_matrix.h5ad"

    if [[ ! -f "$markers_path" || ! -f "$k_path" ]]; then
        echo "Missing markers or K for $sample; generating from ${h5ad_path}"
        python "$PROCESS_DATA" "$sample_dir" --output-dir "$DATA_DIR"
    fi

    k="$(tr -d '[:space:]' < "$k_path")"
    expected_k="$(expected_k_for_sample "$sample")"
    if [[ "$k" != "$expected_k" ]]; then
        echo "Unexpected k for $sample: found $k, expected $expected_k"
        exit 1
    fi

    echo "Launching DLPFC sample $sample with k=$k"

    bash "$REPO_ROOT/experiments/benchmarking/run_benchmark.sh" \
        --data_path="$counts_path" \
        --markers_path="$markers_path" \
        --output_path="$OUTPUT_DIR/$sample" \
        --visium=true \
        --k="$k" \
        --starfysh_lr="1e-6" \
        --hungarian=false &

    sleep 300
done

wait

for sample_dir in "${sample_dirs[@]}"; do
    [ -d "$sample_dir" ] || continue
    sample="$(basename "$sample_dir")"
    h5ad_path="${sample_dir}/${sample}_filtered_feature_bc_matrix.h5ad"
    echo "Analyzing DLPFC results for $sample"
    python "$SCRIPT_DIR/annotate.py" "$h5ad_path" "$OUTPUT_DIR/$sample"
done

python "$SCRIPT_DIR/experiment_results.py" "$OUTPUT_DIR"

echo "DLPFC analysis completed."
