nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/ablation_study/loss_function/run.sh \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/counts_hvgs5000.csv \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/ablation_study/loss_function/TNBC_hvgs5000 \
    5 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/proportions.csv \
    &

nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/ablation_study/loss_function/run.sh \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/mixture_file.csv \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/ablation_study/loss_function/PDAC \
    20 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/proportions.csv \
    &

nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/ablation_study/loss_function/run.sh \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/sdata_lung_s1/pseudospots/pseudospots.csv \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/ablation_study/loss_function/Xenium \
    7 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/sdata_lung_s1/pseudospots/proportions.csv \
    &