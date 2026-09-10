# Gene-Count Ablation

Subsamples genes and runs the full benchmark pipeline at each gene count. Requires SLURM, Python, and all benchmark method environments.

```bash
bash experiments/ablation/gene_count/run.sh \
  data/TNBC/final/TNBC_counts.csv \
  data/TNBC/final/TNBC_marker_genes.csv \
  0.8 \
  experiments/ablation/gene_count/outputs/TNBC \
  5 \
  data/TNBC/final/TNBC_proportions.csv
```

Intermediate subsampled matrices are written to `<output_dir>/tmp` and removed after completion.
