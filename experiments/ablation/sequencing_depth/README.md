# Sequencing-Depth Ablation

Subsamples counts from `0.1` to `1.0` and runs SNMF at each depth. Requires SLURM and the SNMF method environment.

```bash
bash experiments/ablation/sequencing_depth/run.sh \
  data/TNBC/final/TNBC_counts_hvgs5000.csv \
  0.8 \
  experiments/ablation/sequencing_depth/outputs/TNBC \
  5 \
  data/TNBC/final/TNBC_proportions.csv
```

```bash
bash experiments/ablation/sequencing_depth/run.sh \
  data/PDAC/final/PDAC_counts.csv \
  0.8 \
  experiments/ablation/sequencing_depth/outputs/PDAC \
  4 \
  data/PDAC/final/PDAC_proportions.csv
```

```bash
bash experiments/ablation/sequencing_depth/run.sh \
  data/HLC/final/HLC_pseudospots.csv \
  0.8 \
  experiments/ablation/sequencing_depth/outputs/HLC \
  7 \
  data/HLC/final/HLC_proportions.csv
```
