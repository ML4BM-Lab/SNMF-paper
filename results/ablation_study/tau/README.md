nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/tau/run.sh \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/counts.csv \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/tau/TNBC \
    5 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/proportions.csv \
    &

nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/tau/run.sh \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/mixture_file.csv \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/tau/PDAC \
    20 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/proportions.csv \
    &

nohup bash /scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/tau/run.sh \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/bone_marrow_s0/pseudospots/pseudospots.csv \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/results/ablation_study/tau/Xenium \
    5 \
    /scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/bone_marrow_s0/pseudospots/proportions.csv \
    &