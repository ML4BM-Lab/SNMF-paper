# Melanoma Biological Validation

This experiment runs SNMF with `k=4` on the melanoma sample and analyzes the resulting matrices in `main.ipynb`.

Generate SNMF outputs:

```bash
bash methods/SNMF/run.sh \
  data/Melanoma/final/ST_mel1_rep2_counts.csv \
  experiments/melanoma/outputs/K4/ \
  0.5 \
  KL_NB \
  4 \
  "" \
  42 \
  full
```

Then run:

```bash
jupyter notebook experiments/melanoma/main.ipynb
```

The notebook reads `data/Melanoma/final/ST_mel1_rep2.h5ad` and `data/Melanoma/final/melanoma_rep1.png`. Plots and SNMF matrices are stored under `experiments/melanoma/outputs/`.
