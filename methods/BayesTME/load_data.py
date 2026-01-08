import sys

from bayestme import data

import pandas as pd
import scanpy as sc
import numpy as np

import os
import pickle


data_path = sys.argv[1]
output_path = sys.argv[2]

counts = pd.read_csv(data_path, index_col=0)
adata = sc.AnnData(
    X=counts.values.T,
    obs=pd.DataFrame(index=counts.columns),  # rows are observations
    var=pd.DataFrame(index=counts.index)     # rows are variables
)
adata.obs[["array_row", "array_col"]] = adata.obs_names.str.split("x").tolist()

count_mat = adata.X
gene_names = adata.var_names
coordinates = adata.obs[["array_row", "array_col"]].values.astype(int)
n_spots, n_genes = adata.shape

rows = coordinates[:,0]
cols = coordinates[:,1]

n = len(rows)
edges = []

# Build a quick lookup from (row, col) → index
spot_to_idx = {(r, c): i for i, (r, c) in enumerate(zip(rows, cols))}

# Iterate through each spot
for i, (r, c) in enumerate(zip(rows, cols)):
    # Check all possible immediate neighbors (4-connected grid)
    neighbors = [
        (r - 1, c),  # up
        (r + 1, c),  # down
        (r, c - 1),  # left
        (r, c + 1),  # right
    ]
    for nb in neighbors:
        if nb in spot_to_idx:  # if neighbor exists in dataset
            edges.append([i, spot_to_idx[nb]])

edges = np.array(edges, dtype=int)


stdata = data.create_anndata_object(counts=count_mat, 
                                    edges=edges,
                                    coordinates=coordinates, 
                                    tissue_mask=np.ones(n_spots).astype(bool),
                                    gene_names=gene_names, 
                                    layout=data.Layout.SQUARE)
stdata = data.SpatialExpressionDataset(stdata)

with open(os.path.join(output_path, "tmp", "obsnames.pickle"), 'wb') as handle:
    pickle.dump(adata.obs_names, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open(os.path.join(output_path, "tmp", "stdata.pickle"), 'wb') as handle:
    pickle.dump(stdata, handle, protocol=pickle.HIGHEST_PROTOCOL)