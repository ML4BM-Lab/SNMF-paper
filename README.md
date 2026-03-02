nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/run.sh \
    --data_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/mixture_file.csv \
    --markers_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/marker_genes/sig_marker_genes.csv \
    --output_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/results/benchmarking/new/PDAC \
    --k=20 \
    --proportions_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/proportions.csv \
    --starfysh_lr=1e-6 \
    --hungarian=true &
    
nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/run.sh \
    --data_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/counts.csv \
    --markers_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/marker_genes.csv \
    --output_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/results/benchmarking/new/TNBC \
    --k=5 \
    --proportions_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/proportions.csv \
    --starfysh_lr=1e-6 \
    --hungarian=true &