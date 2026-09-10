# PDAC Preparation

Final PDAC experiment inputs live in `data/PDAC/final/` and use:

```text
PDAC_counts.csv
PDAC_marker_genes.csv
PDAC_marker_genes_full.csv
PDAC_proportions.csv
```

Regenerate synthetic counts and proportions into `data/PDAC/final/` with:

```bash
Rscript data/PDAC/scripts/synthetic_data_generation.R
```

Regenerate marker genes into `data/PDAC/final/` with:

```bash
Rscript data/PDAC/scripts/marker_genes.R
```

Generate lower-resolution reannotations with:

```bash
Rscript data/PDAC/scripts/reannotation.R
```

Reannotation scripts read final files from `data/PDAC/final/` and write outputs to `data/PDAC/reannotation/k*/`.
