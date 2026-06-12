nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/ablation_study/ngenes/run.sh \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/counts_hvgs5000.csv \
    0.8 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/ablation_study/ngenes/TNBC \
    5 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/proportions.csv \
    &

nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/ablation_study/ngenes/run.sh \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/mixture_file.csv \
    0.8 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/ablation_study/ngenes/PDAC \
    20 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/proportions.csv \
    &

nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/ablation_study/ngenes/run.sh \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/sdata_lung_s1/pseudospots/pseudospots.csv \
    0.9 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/ablation_study/ngenes/Xenium \
    7 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/sdata_lung_s1/pseudospots/proportions.csv \
    &