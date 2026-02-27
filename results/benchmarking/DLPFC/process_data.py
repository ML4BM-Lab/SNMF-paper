
import sys
import os

import scanpy as sc
import pandas as pd

import tqdm

data_path = sys.argv[1]

adata = sc.read_h5ad(os.path.join(data_path, "filtered_feature_bc_matrix.h5ad"))
# Filter very sparse genes
n_genes = adata.n_vars
sc.pp.filter_genes(adata, min_cells=int(adata.n_obs * 0.1))
print(f"Filtered out {n_genes - adata.n_vars} genes. {adata.n_vars} are kept")

print("Saving counts...")
pd.DataFrame(adata.X.todense().T, index=adata.var_names, columns=[f"{spot.obs['array_row'].values[0]}x{spot.obs['array_col'].values[0]}" for spot in adata]).to_csv(os.path.join(data_path, "counts.csv"))
print("Counts saved correctly!")

if "Region" not in adata.obs.columns:
    print("Assigning cluster ground truth...")
    with open(os.path.join(data_path, f"{os.path.basename(data_path)}_truth.txt"), 'r') as f:
        for line in tqdm.tqdm(f):
            if len(line.strip()) == 2:
                spot, cluster = line.strip().split("\t")
            else:
                continue
            if spot in adata.obs.index:
                adata.obs.loc[spot,"Region"] = cluster
    print("Ground-truth clusters assigned correctly!")

K = adata.obs["Region"].nunique()
with open(os.path.join(data_path, "K.txt"), "w") as f:
    f.write(str(K))
print(f"{K} different clusters in data")

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.tl.rank_genes_groups(
    adata,
    groupby="Region",
    method="wilcoxon"
)

markers = sc.get.rank_genes_groups_df(adata, None)
th = 0.05

clusters_w_sig_markers = -1
while(clusters_w_sig_markers < K):
    if th >= 1:
        print("WARNING: Not all clusters present in marker genes")
        break
    
    sig = markers[
        (markers["pvals_adj"] < th) &
        (markers["logfoldchanges"] > 0)
    ]

    result = sig[["group", "names"]].copy()
    result.columns = ["cluster", "gene"]

    clusters_w_sig_markers = len(list(pd.unique(result["cluster"])))
    th += 0.01

print(f"{clusters_w_sig_markers} different clusters in marker genes")
    
result.to_csv(os.path.join(data_path, "marker_genes.csv"), index=False)