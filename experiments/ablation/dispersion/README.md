# Dispersion Ablation

Compares KL-NB dispersion modes: `full`, `gene`, and `spot`. Requires SLURM and the SNMF method environment.

```bash
bash experiments/ablation/dispersion/run.sh \
  data/TNBC/final/TNBC_counts_hvgs5000.csv \
  experiments/ablation/dispersion/outputs/TNBC \
  5 \
  data/TNBC/final/TNBC_proportions.csv
```

```bash
bash experiments/ablation/dispersion/run.sh \
  data/PDAC/final/PDAC_counts.csv \
  experiments/ablation/dispersion/outputs/PDAC \
  4 \
  data/PDAC/final/PDAC_proportions.csv
```

```bash
bash experiments/ablation/dispersion/run.sh \
  data/HLC/final/HLC_pseudospots.csv \
  experiments/ablation/dispersion/outputs/HLC \
  7 \
  data/HLC/final/HLC_proportions.csv
```
