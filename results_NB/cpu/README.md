nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/cpu/run.sh \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/counts_hvgs5000.csv \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/cpu/results/TNBC \
    5 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/proportions.csv \
    & > /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/cpu/results/TNBC.log

nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/cpu/run.sh \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/mixture_file.csv \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/cpu/results/PDAC \
    20 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/proportions.csv \
    & > /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/cpu/results/PDAC.log

nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/cpu/run.sh \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/sdata_lung_s1/pseudospots/pseudospots.csv \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/cpu/results/Xenium \
    7 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/sdata_lung_s1/pseudospots/proportions.csv \
    & > /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/cpu/results/Xenium.log