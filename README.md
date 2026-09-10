# SNMF Paper Reproducibility Repository

This repository contains the code used to reproduce the experiments for:

> **SNMF: Ultrafast, Spatially-Aware Deconvolution for Spatial Transcriptomics**

Preprint: <https://www.biorxiv.org/content/10.64898/2026.03.17.712043v1.full.pdf>

The `SNMF/` submodule contains the R package implementation. The remaining folders contain data-preparation scripts, method wrappers, and experiment runners for the manuscript.

## Repository Layout

- `SNMF/`: SNMF R package submodule.
- `methods/`: wrappers for SNMF and benchmarked methods.
- `experiments/`: runnable experiment workflows and README files.
- `data/`: data preparation scripts. 
- `assets/`: manuscript figures used in documentation.

Large data files, generated plots, logs, and experiment outputs are ignored by Git. Download the Zenodo archive into each dataset's `final/` folder. DLPFC final files can also be generated from raw per-sample folders with `data/DLPFC/process_data.py`.

## Data

All public data for the manuscript are available from Zenodo:

<https://doi.org/10.5281/zenodo.18852117>

Expected local paths are documented in [data/README.md](data/README.md). The main benchmark inputs are:

- `data/TNBC/final/TNBC_counts_hvgs5000.csv`
- `data/PDAC/final/PDAC_counts.csv`
- `data/HLC/final/HLC_pseudospots.csv`
- `data/Melanoma/final/ST_mel1_rep2_counts.csv`
- `data/DLPFC/final/{sample_id}/{sample_id}_counts.csv`

## Environment

The manuscript experiments were run on a SLURM-managed HPC cluster with **R 4.4.1**, **Python 3.9**, Apptainer, and NVIDIA RTX 3090 GPU nodes for GPU benchmarks. Install Python dependencies with:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Install R dependencies as needed for the selected methods:

```r
install.packages("devtools")
devtools::install_github("YMa-lab/CARD")
devtools::install_github("qunhualilab/retrofit")
devtools::install_github("yyolanda/SMART")
install.packages("BiocManager")
BiocManager::install("STdeconvolve")
install.packages("RcppHungarian")
```

## Reproducing Experiments

### Benchmarks

Benchmark TNBC (**Figure 1A-E** in the paper):

```bash
bash experiments/benchmarking/run_benchmark.sh \
  --data_path=data/TNBC/final/TNBC_counts_hvgs5000.csv \
  --markers_path=data/TNBC/final/TNBC_marker_genes.csv \
  --output_path=experiments/benchmarking/outputs/TNBC \
  --k=5 \
  --proportions_path=data/TNBC/final/TNBC_proportions.csv \
  --hungarian=true
```

Benchmark all DLPFC samples (**Figure 1F** and **Supp. Figure 11** in the paper):

```bash
bash experiments/benchmarking/dlpfc/run.sh
```

Benchmark PDAC (**Supp. Figures 6 and 8-9** in the paper):

```bash
bash experiments/benchmarking/run_benchmark.sh \
  --data_path=data/TNBC/final/TNBC_counts_hvgs5000.csv \
  --markers_path=data/TNBC/final/TNBC_marker_genes.csv \
  --output_path=experiments/benchmarking/outputs/TNBC \
  --k=5 \
  --proportions_path=data/TNBC/final/TNBC_proportions.csv \
  --hungarian=true
```

Benchmark HLC (**Supp. Figures 7 and 10** in the paper):

```bash
bash experiments/benchmarking/run_benchmark.sh \
  --data_path=data/TNBC/final/TNBC_counts_hvgs5000.csv \
  --markers_path=data/TNBC/final/TNBC_marker_genes.csv \
  --output_path=experiments/benchmarking/outputs/TNBC \
  --k=5 \
  --proportions_path=data/TNBC/final/TNBC_proportions.csv \
  --hungarian=true
```

### Ablation studies

Effect of the spatial mixing matrix $\mathbf{S}$ on deconvolution accuracy: performance across different values of $\tau$ (**Supp. Figure 1** in the paper):

```bash
bash experiments/ablation/tau/run.sh \
  data/TNBC/final/TNBC_counts_hvgs5000.csv \
  experiments/ablation/tau/outputs/TNBC \
  5 \
  data/TNBC/final/TNBC_proportions.csv

bash experiments/ablation/tau/run.sh \
  data/PDAC/final/PDAC_counts.csv \
  experiments/ablation/tau/outputs/PDAC \
  20 \
  data/PDAC/final/PDAC_proportions.csv

bash experiments/ablation/tau/run.sh \
  data/HLC/final/HLC_counts.csv \
  experiments/ablation/tau/outputs/HLC \
  7 \
  data/HLC/final/HLC_proportions.csv
```

Effect of the loss function to optimize $\mathcal{L}$ on deconvolution accuracy: performance for spot-wise *Frobenius norm*, *Poisson Kullback-Leibler divergence* and *Negative Binomial negative log-likelihood* as optimization alternatives (**Supp. Figure 2** in the paper):

```bash
bash experiments/ablation/loss_function/run.sh \
  data/TNBC/final/TNBC_counts_hvgs5000.csv \
  experiments/ablation/loss_function/outputs/TNBC \
  5 \
  data/TNBC/final/TNBC_proportions.csv

bash experiments/ablation/loss_function/run.sh \
  data/PDAC/final/PDAC_counts.csv \
  experiments/ablation/loss_function/outputs/PDAC \
  20 \
  data/PDAC/final/PDAC_proportions.csv

bash experiments/ablation/loss_function/run.sh \
  data/HLC/final/HLC_counts.csv \
  experiments/ablation/loss_function/outputs/HLC \
  7 \
  data/HLC/final/HLC_proportions.csv
```

Effect of the number of *highly variable genes (HVGs)* on deconvolution accuracy for the TNBC dataset (**Supp. Figure 3** in the paper):

```bash
bash experiments/ablation/gene_count/run.sh \
  data/TNBC/final/TNBC_counts_hvgs5000.csv \
  experiments/ablation/gene_count/outputs/TNBC \
  5 \
  data/TNBC/final/TNBC_proportions.csv

bash experiments/ablation/gene_count/run.sh \
  data/PDAC/final/PDAC_counts.csv \
  experiments/ablation/gene_count/outputs/PDAC \
  20 \
  data/PDAC/final/PDAC_proportions.csv

bash experiments/ablation/gene_count/run.sh \
  data/HLC/final/HLC_counts.csv \
  experiments/ablation/gene_count/outputs/HLC \
  7 \
  data/HLC/final/HLC_proportions.csv
```

### Runtime comparison (Supp. Figure  in the paper)

```bash
bash experiments/runtime/run.sh \
  data/TNBC/final/TNBC_counts_hvgs5000.csv \
  experiments/runtime/outputs/results/TNBC \
  5 \
  data/TNBC/final/TNBC_proportions.csv
```

### Melanoma biological validation (Figure 2 and Supp. Figures 12-19 in the paper)

```bash
bash methods/SNMF/run.sh \
  data/Melanoma/final/ST_mel1_rep2_counts.csv \
  experiments/melanoma/outputs/K4/ \
  0.5 \
  KL_NB \
  4 \
  "" \
  42 \
  full
```

Then run [experiments/melanoma/main.ipynb](experiments/melanoma/main.ipynb).
