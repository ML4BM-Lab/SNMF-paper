nohup bash -c '
    cd /scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SNMF
    mkdir -p "/scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/melanoma/K4"
    bash run.sh \
        "/scratch/lalonsoeste/PhD/SpatialTranscriptomics/data/spatial/BayesSpace/melanoma_thrane/processed/ST_mel1_rep2_counts.csv" \
        "/scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/melanoma/K4/" \
        "0.8" \
        "KL_NB" \
        "4" \
        "" \
        "42"
' > "/scratch/lalonsoeste/PhD/NMF_deconvolution/results_NB/melanoma/K4.log" 2>&1 &