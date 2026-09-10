# Benchmarking

Runs SNMF, NMF, CARD, RETROFIT, STdeconvolve, SMART, Starfysh, BayesTME, and SpiceMix on one dataset.

The workflow requires SLURM, the cluster R/Python modules used by `methods/*`, and a populated `data/` directory. Outputs are written to `experiments/benchmarking/outputs/`.

## Command Template

```bash
bash experiments/benchmarking/run_benchmark.sh \
  --data_path=<counts.csv> \
  --markers_path=<marker_genes.csv> \
  --output_path=<output_dir> \
  --k=<components> \
  --proportions_path=<ground_truth_proportions.csv> \
  --hungarian=true
```

## Manuscript Datasets

```bash
bash experiments/benchmarking/run_benchmark.sh \
  --data_path=data/TNBC/final/TNBC_counts_hvgs5000.csv \
  --markers_path=data/TNBC/final/TNBC_marker_genes.csv \
  --output_path=experiments/benchmarking/outputs/TNBC \
  --k=5 \
  --proportions_path=data/TNBC/final/TNBC_proportions.csv \
  --hungarian=true
```

```bash
bash experiments/benchmarking/run_benchmark.sh \
  --data_path=data/PDAC/final/PDAC_counts.csv \
  --markers_path=data/PDAC/final/PDAC_marker_genes.csv \
  --output_path=experiments/benchmarking/outputs/PDAC \
  --k=4 \
  --proportions_path=data/PDAC/final/PDAC_proportions.csv \
  --hungarian=true
```

```bash
bash experiments/benchmarking/run_benchmark.sh \
  --data_path=data/HLC/final/HLC_pseudospots.csv \
  --markers_path=data/HLC/final/HLC_marker_genes.csv \
  --output_path=experiments/benchmarking/outputs/HLC \
  --k=7 \
  --proportions_path=data/HLC/final/HLC_proportions.csv \
  --hungarian=true
```

## Friedman-Nemenyi Analysis

Run this after the benchmark outputs exist:

```bash
python3 experiments/benchmarking/scripts/friedman_nemenyi.py experiments/benchmarking/outputs
```
