# Loss-Function Ablation

Runs SNMF with squared-error, KL-Poisson, and KL-NB losses. Requires SLURM and the SNMF method environment.

```bash
bash experiments/ablation/loss_function/run.sh \
  data/TNBC/final/TNBC_counts_hvgs5000.csv \
  experiments/ablation/loss_function/outputs/TNBC \
  5 \
  data/TNBC/final/TNBC_proportions.csv
```

```bash
bash experiments/ablation/loss_function/run.sh \
  data/PDAC/final/PDAC_counts.csv \
  experiments/ablation/loss_function/outputs/PDAC \
  4 \
  data/PDAC/final/PDAC_proportions.csv
```

```bash
bash experiments/ablation/loss_function/run.sh \
  data/HLC/final/HLC_pseudospots.csv \
  experiments/ablation/loss_function/outputs/HLC \
  7 \
  data/HLC/final/HLC_proportions.csv
```
