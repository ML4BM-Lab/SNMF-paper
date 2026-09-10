# TNBC Preparation

Final TNBC experiment inputs live in `data/TNBC/final/` and use:

```text
TNBC_counts.csv
TNBC_counts_hvgs5000.csv
TNBC_marker_genes.csv
TNBC_proportions.csv
```

Generate the HVG-filtered matrix with:

```bash
Rscript data/TNBC/scripts/hvgs.R
```

Generate spatial subsets for runtime experiments with:

```bash
python3 data/TNBC/scripts/subset.py
```

Generated subset matrices are written to `data/TNBC/subsets/`.
