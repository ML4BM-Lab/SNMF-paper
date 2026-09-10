# DLPFC Data

DLPFC final experiment inputs live under `data/DLPFC/final/`, with one folder per sample.

```text
data/DLPFC/final/{sample_id}/
```

Each sample folder contains the Zenodo-style `{sample_id}_` files used by `experiments/benchmarking/dlpfc/run.sh`.

Generate the final folder from local per-sample inputs with:

```bash
python3 data/DLPFC/process_data.py /path/to/DLPFC
```

The command accepts either a raw DLPFC root containing sample folders or a single raw sample folder.
