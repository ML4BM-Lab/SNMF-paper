# Tau Ablation

Sweeps SNMF spatial regularization values from `0.1` to `1.0`. Requires SLURM and the SNMF method environment.

```bash
bash experiments/ablation/tau/run.sh \
  data/TNBC/final/TNBC_counts_hvgs5000.csv \
  experiments/ablation/tau/outputs/TNBC \
  5 \
  data/TNBC/final/TNBC_proportions.csv
```

```bash
bash experiments/ablation/tau/run.sh \
  data/PDAC/final/PDAC_counts.csv \
  experiments/ablation/tau/outputs/PDAC \
  4 \
  data/PDAC/final/PDAC_proportions.csv
```

```bash
bash experiments/ablation/tau/run.sh \
  data/HLC/final/HLC_pseudospots.csv \
  experiments/ablation/tau/outputs/HLC \
  7 \
  data/HLC/final/HLC_proportions.csv
```

Compare tau signatures against marker genes after outputs exist:

```bash
python3 experiments/ablation/tau/signatures.py --results-path experiments/ablation/tau/outputs --dataset all
```
