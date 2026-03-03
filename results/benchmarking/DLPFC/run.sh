#!/bin/bash

base_dir="/scratch/lalonsoeste/PhD/SpatialTranscriptomics/data/spatial/10x/Visium/DLPFC"

for folder in "$base_dir"/*; do
    if [ -d "$folder" ]; then
        sample=$(basename "$folder")

        counts_path="${folder}/counts.csv"
        markers_path="${folder}/marker_genes.csv"
        K_path="${folder}/K.txt"

        echo "🔎 Processing sample: $sample"

        # Step 1 — Generate data if missing
        if [[ ! -f "$counts_path" || ! -f "$markers_path" || ! -f "$K_path" ]]; then
            echo "⚙️  Missing counts or markers for $sample"
            echo "   → Running process_data.py"

            source /scratch/lalonsoeste/PhD/NMF_deconvolution/.venv/bin/activate
            python /scratch/lalonsoeste/PhD/NMF_deconvolution/results/benchmarking/DLPFC/process_data.py "$folder"

            echo "✅ Data generated for $sample"
        fi

        # Step 2 — Launch pipeline (always after generation check)
        echo "🚀 Launching pipeline for $sample"

        k=$(cat "$folder/K.txt")

        bash /scratch/lalonsoeste/PhD/NMF_deconvolution/run.sh \
            --data_path="$counts_path" \
            --markers_path="$markers_path" \
            --output_path="/scratch/lalonsoeste/PhD/NMF_deconvolution/results/benchmarking/DLPFC/new/${sample}" \
            --visium=true \
            --k="$k" \
            --starfysh_lr="1e-6" \
            --hungarian=false &

        # Wait 5 minutes before next launch
        sleep 300

    fi
done

wait

echo "🎉 All samples processed."

source /scratch/lalonsoeste/PhD/NMF_deconvolution/.venv/bin/activate
for folder in "$base_dir"/*; do
    if [ -d "$folder" ]; then
        sample=$(basename "$folder")

        echo "🚀 Analyzing results for $sample"

        python /scratch/lalonsoeste/PhD/NMF_deconvolution/results/benchmarking/DLPFC/annotate.py \
            "${folder}/filtered_feature_bc_matrix.h5ad" \
            "/scratch/lalonsoeste/PhD/NMF_deconvolution/results/benchmarking/DLPFC/new/${sample}"
    fi
done

echo "🎉 Analysis completed."
