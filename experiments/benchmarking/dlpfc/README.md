# DLPFC Benchmarking

Runs the benchmark pipeline across all DLPFC Visium samples in `data/DLPFC/final`.

Each sample has its own folder:

```text
data/DLPFC/final/{sample_id}/
```

Each sample folder is expected to have `{sample_id}_counts.csv`, `{sample_id}_marker_genes.csv`, `{sample_id}_K.txt`, `{sample_id}_truth.txt`, and `{sample_id}_filtered_feature_bc_matrix.h5ad`. The runner reads `{sample_id}_K.txt` to set `k`.

To generate the final folders from raw per-sample DLPFC inputs, run:

```bash
python3 data/DLPFC/process_data.py /path/to/DLPFC
```

```bash
bash experiments/benchmarking/dlpfc/run.sh
```

Per-sample method outputs are written to `experiments/benchmarking/outputs/DLPFC/{sample_id}`. Aggregate ARI plots and tables are written under `experiments/benchmarking/outputs/DLPFC`.
