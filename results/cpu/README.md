nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results/cpu/run.sh \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/counts.csv \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/results/cpu/results/TNBC \
    5 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/proportions.csv \
    & > /scratch/lalonsoeste/PhD/NMF_deconvolution/results/cpu/results/TNBC.log

nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results/cpu/run.sh \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/mixture_file.csv \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/results/cpu/results/PDAC \
    20 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/proportions.csv \
    & > /scratch/lalonsoeste/PhD/NMF_deconvolution/results/cpu/results/PDAC.log

nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results/cpu/run.sh \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/sdata_lung_s1/pseudospots/pseudospots.csv \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/results/cpu/results/Xenium \
    7 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/sdata_lung_s1/pseudospots/proportions.csv \
    & > /scratch/lalonsoeste/PhD/NMF_deconvolution/results/cpu/results/Xenium.log