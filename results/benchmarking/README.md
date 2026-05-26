nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results/benchmarking/benchmark.sh \
    --data_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/counts.csv \
    --markers_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/marker_genes.csv \
    --output_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/results/benchmarking/TNBC \
    --k=5 \
    --proportions_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/proportions.csv \
    --starfysh_lr=1e-6 \
    & > /scratch/lalonsoeste/PhD/NMF_deconvolution/results/benchmarking/TNBC.log

nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results/benchmarking/benchmark.sh \
    --data_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/mixture_file.csv \
    --markers_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/marker_genes/marker_genes.csv \
    --output_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/results/benchmarking/PDAC \
    --k=20 \
    --proportions_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/proportions.csv \
    --starfysh_lr=1e-6 \
    & > /scratch/lalonsoeste/PhD/NMF_deconvolution/results/benchmarking/PDAC.log

nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results/benchmarking/benchmark.sh \
    --data_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/bone_marrow_s0/pseudospots/pseudospots.csv \
    --markers_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/bone_marrow_s0/pseudospots/marker_genes.csv \
    --visium=true \
    --output_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/results/benchmarking/Xenium \
    --k=5 \
    --proportions_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/bone_marrow_s0/pseudospots/proportions.csv \
    --starfysh_lr=1e-6 \
    & > /scratch/lalonsoeste/PhD/NMF_deconvolution/results/benchmarking/Xenium.log