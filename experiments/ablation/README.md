# Ablation Studies

Each subdirectory contains one SNMF ablation workflow. Outputs are written under that experiment's `outputs/` directory and ignored by Git.

Run commands from the repository root. The wrappers normalize relative paths internally.

- `loss_function/`: compares squared error, KL-Poisson, and KL-NB losses.
- `tau/`: sweeps the spatial regularization parameter.
- `dispersion/`: compares KL-NB dispersion modes.
- `gene_count/`: benchmarks performance across gene subset sizes.
- `sequencing_depth/`: benchmarks performance across read-depth subsampling levels.
