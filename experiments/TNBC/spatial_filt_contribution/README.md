# scDesign3
nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/experiments/starfysh_synthetic/spatial_filt_contribution/run.sh \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/mixture_file.csv \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/experiments/starfysh_synthetic/spatial_filt_contribution/results/scDesign3/ \
    20 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/proportions.csv \
    > /scratch/lalonsoeste/PhD/NMF_deconvolution/experiments/starfysh_synthetic/spatial_filt_contribution/scDesign3.out &

# Starfysh
nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/experiments/starfysh_synthetic/spatial_filt_contribution/run.sh \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/starfysh/simulated/processed/counts.csv \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/experiments/starfysh_synthetic/spatial_filt_contribution/results/starfysh/ \
    5 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/starfysh/simulated/processed/proportions.csv \
    > /scratch/lalonsoeste/PhD/NMF_deconvolution/experiments/starfysh_synthetic/spatial_filt_contribution/starfysh.out &

