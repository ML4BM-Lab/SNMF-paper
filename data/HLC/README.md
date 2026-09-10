# HLC Preparation

HLC contains pseudospots generated from the STHELAR human lung cancer data. Final experiment inputs live in `data/HLC/final/` and use:

```text
HLC_pseudospots.csv
HLC_marker_genes.csv
HLC_marker_genes_full.csv
HLC_proportions.csv
```

Download the source zarr archive with:

```bash
bash data/HLC/scripts/download.sh
```

Generate pseudospots with:

```bash
cd data/HLC/scripts
sbatch pseudospots_generation.slurm
```

Final CSV files are written to `data/HLC/final/`. Derived h5ad files, plots, and SLURM output logs are written to `data/HLC/derived/`.
