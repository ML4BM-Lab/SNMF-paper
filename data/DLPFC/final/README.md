# DLPFC Data

DLPFC samples use the `{sample_id}_` prefix. Each sample has its own folder:

```text
data/DLPFC/final/{sample_id}/
```

Each sample folder must include:

```text
{sample_id}_counts.csv
{sample_id}_marker_genes.csv
{sample_id}_K.txt
{sample_id}_truth.txt
{sample_id}_filtered_feature_bc_matrix.h5ad
```

`{sample_id}_K.txt` contains one integer. The DLPFC benchmark reads it to set `k` automatically.

Generate this directory from local per-sample folders with:

```bash
python3 data/DLPFC/process_data.py /path/to/DLPFC
```
