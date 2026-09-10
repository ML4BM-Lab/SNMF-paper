# Runtime Comparison

Compares SNMF CPU and GPU runtimes across seeds. This workflow requires SLURM, Apptainer, and the SNMF container image expected by `methods/SNMF`.

Outputs are written under `experiments/runtime/outputs/`.

```bash
bash experiments/runtime/run.sh \
  data/TNBC/final/TNBC_counts_hvgs5000.csv \
  experiments/runtime/outputs/TNBC \
  5 \
  data/TNBC/final/TNBC_proportions.csv
```

```bash
bash experiments/runtime/run.sh \
  data/PDAC/final/PDAC_counts.csv \
  experiments/runtime/outputs/PDAC \
  4 \
  data/PDAC/final/PDAC_proportions.csv
```

```bash
bash experiments/runtime/run.sh \
  data/HLC/final/HLC_pseudospots.csv \
  experiments/runtime/outputs/HLC \
  7 \
  data/HLC/final/HLC_proportions.csv
```
