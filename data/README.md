# Data

This directory separates final experiment inputs from dataset-preparation material inside each dataset.

- `data/{dataset}/final/`: canonical local mirror of the final Zenodo-style files used by experiments.
- `data/{dataset}/scripts/`: dataset-preparation scripts when applicable.
- Other dataset-specific folders contain source files and generated intermediate outputs.

Large data files are ignored by Git. Populate each `final/` folder from the [Zenodo](<https://doi.org/10.5281/zenodo.18852117>) archive before running experiments

Zenodo filename prefixes are preserved in each final dataset folder:

- `PDAC_`: pancreatic ductal adenocarcinoma synthetic spatial mixture.
- `TNBC_`: triple-negative breast cancer synthetic spatial mixture.
- `ST_mel1_rep2`: melanoma validation dataset.
- `{sample_id}_`: DLPFC Visium samples, for example `151673_counts.csv`.

For DLPFC, `{sample_id}_K.txt` contains the number of annotated regions and is read automatically by `experiments/benchmarking/dlpfc/run.sh`.

## Expected Files

```text
data/
  TNBC/
    final/
      TNBC_counts.csv
      TNBC_counts_hvgs5000.csv
      TNBC_marker_genes.csv
      TNBC_proportions.csv
  PDAC/
    final/
      PDAC_counts.csv
      PDAC_marker_genes.csv
      PDAC_marker_genes_full.csv
      PDAC_proportions.csv
  HLC/
    final/
      HLC_pseudospots.csv
      HLC_marker_genes.csv
      HLC_marker_genes_full.csv
      HLC_proportions.csv
  Melanoma/
    final/
      ST_mel1_rep2_counts.csv
      ST_mel1_rep2.h5ad
      melanoma_rep1.png
  DLPFC/
    final/
      {sample_id}/
        {sample_id}_counts.csv
        {sample_id}_marker_genes.csv
        {sample_id}_K.txt
        {sample_id}_truth.txt
        {sample_id}_filtered_feature_bc_matrix.h5ad
```

## Generate DLPFC Final Files

```bash
python3 data/DLPFC/process_data.py /path/to/DLPFC
```

The DLPFC script accepts either a root containing per-sample folders or one sample folder. It writes final files into `data/DLPFC/final/{sample_id}/`.
