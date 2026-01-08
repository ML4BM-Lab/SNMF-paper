
import sys
import os

import pandas as pd
import numpy as np
import scanpy as sc

from starfysh import (AA, utils)

import pickle


data_path = sys.argv[1]
marker_genes_path = sys.argv[2]
output_path = sys.argv[3]

## Load data
data = pd.read_csv(data_path, index_col=0)
adata_raw = sc.AnnData(
    X=data.values.T,
    obs=pd.DataFrame(index=data.columns),  # rows are observations
    var=pd.DataFrame(index=data.index)     # rows are variables
)
adata_raw.obs[["array_row", "array_col"]] = adata_raw.obs_names.str.split("x").tolist()

# Filter
adata = adata_raw.copy()

adata.var['mt'] = np.logical_or(
    adata.var_names.str.startswith('MT-'),
    adata.var_names.str.startswith('mt-')
)
adata.var['rb'] = adata.var_names.str.startswith(('RP', 'Rp', 'rp'))

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
mask_cell = adata.obs['pct_counts_mt'] < 100
mask_gene = np.logical_and(~adata.var['mt'], ~adata.var['rb'])

adata = adata[mask_cell, mask_gene]
sc.pp.filter_genes(adata, min_cells=1)

sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=2000, inplace=True)

adata_raw = adata_raw[adata.obs_names, adata.var_names]
adata_raw.var['highly_variable'] = adata.var['highly_variable']
adata_raw.obs = adata.obs

# Marker genes
gene_sig = pd.read_csv(marker_genes_path)
markers = gene_sig.groupby("cluster")["gene"].apply(list)
max_len = markers.apply(len).max()
gene_sig = pd.DataFrame({ct: genes + [None]*(max_len - len(genes))
                        for ct, genes in markers.items()})
gene_sig = utils.filter_gene_sig(gene_sig, adata.to_df())

# Image metadata
map_info = adata.obs[["array_row", "array_col"]]
img_metadata = {
        'img': None,
        'map_info': map_info,
        'scalefactor': None
    }

# Finding anchor spots
visium_args = utils.VisiumArguments(adata_raw,
                                    adata,
                                    gene_sig,
                                    img_metadata,
                                    window_size=3, # adjust window_size for considering the neighbor density
                                    sample_id="sample_id")

## Archetypal analysis
aa_model = AA.ArchetypalAnalysis(adata_orig=adata)
try:
    archetype, arche_dict, major_idx, evs = aa_model.compute_archetypes()
except:
    archetype, arche_dict, major_idx, evs = aa_model.compute_archetypes(n_iters=200)

# (1). Find archetypal spots & archetypal clusters
arche_df = aa_model.find_archetypal_spots(major=True)

# (2). Find marker genes associated with each archetypal cluster
markers_df = aa_model.find_markers(display=False)

# (3). Map archetypes to the closest anchors within `r` nearest neighbors
# Choose the top `anchor_percent` (N%) anchors per cell type for archetype mapping
# In general, set lower `anchor_percent` for fine resolved cell-states
anchors_df = visium_args.get_anchors()
anchor_percent = 0.05
n_top_anchors = int(anchor_percent*adata_raw.shape[0])
map_df, map_dict = aa_model.assign_archetypes(anchor_df=anchors_df[:n_top_anchors],
                                              r=n_top_anchors)

visium_args = utils.refine_anchors(visium_args,
                                   aa_model,
                                   thld=anchor_percent,
                                   n_genes=10)

# Write data
adata_raw.write_h5ad(os.path.join(output_path, "tmp", "adata.h5ad"))
with open(os.path.join(output_path, "tmp", "visium_args.pickle"), 'wb') as handle:
    pickle.dump(visium_args, handle, protocol=pickle.HIGHEST_PROTOCOL)
gene_sig.to_csv(os.path.join(output_path, "tmp", "gene_sig.csv"))