#!/bin/bash

base_dir="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/starfysh/simulated/subsets"

cd "/scratch/lalonsoeste/PhD/NMF_deconvolution/experiments/starfysh_synthetic/spots_subsets"

for folder in "$base_dir"/*; do
    if [ -d "$folder" ]; then
        folder_name=$(basename "$folder")

        echo "Pipeline launched for $folder_name"

        bash /scratch/lalonsoeste/PhD/NMF_deconvolution/run.sh \
            --data_path="$folder/counts.csv" \
            --markers_path="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/starfysh/simulated/processed/marker_genes.csv" \
            --output_path="/scratch/lalonsoeste/PhD/NMF_deconvolution/experiments/starfysh_synthetic/spots_subsets/results/$folder_name" \
            --k=5 \
            --proportions_path="$folder/proportions.csv" \
            --starfysh_lr="1e-6" \
            --hungarian=true &

        sleep 300
    fi
done

wait  # final wait for any leftover jobs

# Plot results
python ./plot.py ./results/
