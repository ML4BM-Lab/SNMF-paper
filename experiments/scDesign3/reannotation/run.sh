#!/bin/bash

base_dir="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/reannotation"

for folder in "$base_dir"/k*/; do
    if [ -d "$folder" ]; then
        # Extract filename without path
        folder_name=$(basename "$folder")

        # Extract k value from folder name (e.g., k5 → 5)
        k=$(echo "$folder_name" | sed -E 's/k([0-9]+)/\1/')

        proportions_path="${folder}/proportions_k${k}.csv"
        markers_path="${folder}/marker_genes_k${k}.csv"

        if [[ -f "$proportions_path" && -f "$markers_path" ]]; then
            echo "🚀 Launching pipeline for $k cell types"

            bash /scratch/lalonsoeste/PhD/NMF_deconvolution/run.sh \
                --data_path="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/mixture_file.csv" \
                --markers_path="$markers_path" \
                --output_path="/scratch/lalonsoeste/PhD/NMF_deconvolution/experiments/scDesign3/reannotation/results/$k" \
                --k="$k" \
                --proportions_path="$proportions_path" \
                --starfysh_lr="1e-6" \
                --hungarian=true &

            # Wait 5 minutes before launching the next job
            sleep 300
        else
            echo "⚠️ Missing files for $folder_name — skipping."
        fi
    fi
done

wait  # Wait for any background jobs to finish

# Plot results afterwards
source /scratch/lalonsoeste/PhD/NMF_deconvolution/.venv/bin/activate
python /scratch/lalonsoeste/PhD/NMF_deconvolution/experiments/scDesign3/reannotation/plot.py \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/experiments/scDesign3/reannotation/results/ \
    $base_dir \
    "true"
