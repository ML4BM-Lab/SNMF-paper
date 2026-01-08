nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/experiments/starfysh_synthetic/spots_subsets/run.sh &

nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/run.sh \
    --data_path="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/starfysh/simulated/subsets/2116spots_48.7x48.7/counts.csv" \
    --markers_path="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/starfysh/simulated/processed/marker_genes.csv" \
    --output_path="/scratch/lalonsoeste/PhD/NMF_deconvolution/experiments/starfysh_synthetic/spots_subsets/new/results/2116spots_48.7x48.7" \
    --k=5 \
    --proportions_path="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/starfysh/simulated/subsets/2116spots_48.7x48.7/proportions.csv" \
    --hungarian=true &