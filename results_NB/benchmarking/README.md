nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/benchmarking/benchmark.sh \
    --data_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/counts_hvgs5000.csv \
    --markers_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/marker_genes.csv \
    --output_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/benchmarking/TNBC \
    --k=5 \
    --proportions_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/proportions.csv \
    --snmf_tau=0.8 \
    --starfysh_lr=1e-6 \
    & > /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/benchmarking/TNBC.log

nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/benchmarking/benchmark.sh \
    --data_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/mixture_file.csv \
    --markers_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/marker_genes/marker_genes.csv \
    --output_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/benchmarking/PDAC \
    --k=20 \
    --proportions_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/proportions.csv \
    --snmf_tau=0.8 \
    --starfysh_lr=1e-6 \
    & > /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/benchmarking/PDAC.log

nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/benchmarking/benchmark.sh \
    --data_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/sdata_lung_s1/pseudospots/pseudospots.csv \
    --markers_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/sdata_lung_s1/pseudospots/marker_genes.csv \
    --visium=true \
    --output_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/benchmarking/Xenium \
    --k=7 \
    --proportions_path=/scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/sdata_lung_s1/pseudospots/proportions.csv \
    --snmf_tau=0.9 \
    --starfysh_lr=1e-6 \
    & > /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/benchmarking/Xenium.log

python /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/benchmarking/friedman_nemenyi/friedman_nemenyi.py /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/benchmarking