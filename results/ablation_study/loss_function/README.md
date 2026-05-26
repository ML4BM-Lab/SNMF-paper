nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/loss_function/run.sh \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/counts.csv \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/loss_function/TNBC \
    5 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/proportions.csv \
    &

nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/loss_function/run.sh \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/mixture_file.csv \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/loss_function/PDAC \
    20 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/proportions.csv \
    &

nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/loss_function/run.sh \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/bone_marrow_s0/pseudospots/pseudospots.csv \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/loss_function/Xenium \
    5 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/bone_marrow_s0/pseudospots/proportions.csv \
    &