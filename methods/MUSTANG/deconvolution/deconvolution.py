#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np

from pathlib import Path

from bayestme import data, gene_filtering
from src import data

data_path = sys.argv[1]
output_path = sys.argv[2]
technology = sys.argv[3]

if ".csv" in data_path:
    counts = pd.read_csv(data_path, index_col=0)
else:
    raise ValueError

raw_counts = counts.values
pos = np.array(counts.columns.str.split("x").tolist())
genes = counts.index.values
mask = [1] * counts.shape[1]
barcodes = counts.columns

similarity_graph = pd.read_csv(f"{output_path}/similarity_neighbors.csv", index_col=0).reset_index(drop=True).to_numpy()

stdata = data.SpatialExpressionDataset.from_arrays(
    raw_counts=raw_counts.T,
    positions=pos,
    tissue_mask=mask,
    gene_names=genes,
    layout=data.Layout.HEX if technology == "Visium" else data.Layout.SQUARE,
    barcodes=barcodes
)

stdata_stddev = gene_filtering.select_top_genes_by_standard_deviation(stdata, n_gene=1000)
stdata_threshold = gene_filtering.filter_genes_by_spot_threshold(stdata_stddev, spot_threshold=0.95)
stdata_filtered = gene_filtering.filter_ribosome_genes(stdata_threshold)

best_lambda = 1000
best_n_components = 11
deconv_res = deconvolution.deconvolve(
    reads=stdata_filtered.reads,
    edges=stdata_filtered.edges,
    n_gene=1000,
    n_components=best_n_components,
    lam2=best_lambda,
    n_samples=100,
    n_burnin=1000,
    n_thin=5,
    bkg=False,
    lda=False,
    similarity_graph=similarity_graph
)

Path(f"{output_path}/deconvolution_plots_transcrip").mkdir(exist_ok=True)
Path(f"{output_path}/deconvolution_res_transcrip").mkdir(exist_ok=True)

deconvolution.add_deconvolution_results_to_dataset(stdata_filtered, deconv_res)
deconv_res.save(f'{output_path}/deconvolution_res_transcrip/deconvolve_result.h5')

deconvolution.plot_cell_num_scatterpie(stdata_filtered, f'{output_path}/deconvolution_plots_transcrip/Scatter_piechart.png')
deconvolution.plot_cell_prob(stdata_filtered, output_dir=f'{output_path}/deconvolution_plots_transcrip', output_format='png')
deconvolution.plot_cell_num(stdata_filtered, output_dir=f'{output_path}/deconvolution_plots_transcrip', output_format='png')
