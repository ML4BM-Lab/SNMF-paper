import sys
import os
import itertools

import scanpy as sc
import pandas as pd

from sklearn.metrics import adjusted_rand_score

data_path = sys.argv[1]
results_path = sys.argv[2]

adata = sc.read_h5ad(data_path)

methods = []

# -------------------------------------------------
# Load predictions
# -------------------------------------------------
for (dirpath, dirnames, filenames) in os.walk(results_path):
    for file in filenames:
        if (
            file.endswith(".csv") and 
            "proportions" in file and 
            "tmp" not in dirpath
        ):
            method = dirpath.split("/")[-1]

            proportions = pd.read_csv(
                os.path.join(dirpath, file),
                index_col=0
            )
            if proportions.shape[0] < proportions.shape[1]:
                proportions = proportions.transpose()

            # Assign major predicted cell type
            adata.obs["rowcol"] = (
                adata.obs["array_row"].astype(str)
                + "x"
                + adata.obs["array_col"].astype(str)
            )

            major_ct = proportions.idxmax(axis="columns")
            mapping = major_ct.to_dict()
            adata.obs[f"{method}_major_ct"] = adata.obs["rowcol"].map(mapping)

            # adata.obs[f"{method}_major_ct"] = adata.obs[f"{method}_major_ct"].astype("category")

            methods.append(method)

methods = list(set(methods))  # remove duplicates

# -------------------------------------------------
# Plot predictions
# -------------------------------------------------
import matplotlib.pyplot as plt

plots_dir = os.path.join(results_path, "plots")
os.makedirs(plots_dir, exist_ok=True)

# Ground truth
adata.obs["Region"] = adata.obs["Region"].astype("category")
sc.pl.spatial(
    adata,
    color="Region",
    show=False,
    title="Ground truth"
)
plt.savefig(os.path.join(plots_dir, "ground_truth.png"),
            dpi=300,
            bbox_inches="tight")
plt.close()

# Methods
for method in methods:
    sc.pl.spatial(
        adata,
        color=f"{method}_major_ct",
        show=False,
        title=method
    )
    plt.savefig(os.path.join(plots_dir, f"{method}.png"),
                dpi=300,
                bbox_inches="tight")
    plt.close()

# -------------------------------------------------
# Compute ARIs
# -------------------------------------------------
aris = {}

true_labels = adata.obs["Region"].astype(str)

for method in methods:

    pred_labels = adata.obs[f"{method}_major_ct"].astype(str)

    # Keep only cells with predictions
    mask = pred_labels.notna()
    y_true = true_labels[mask]
    y_pred = pred_labels[mask]

    # Case 1: labels already match reference vocabulary
    if set(y_pred.unique()).issubset(set(y_true.unique())):
        ari = adjusted_rand_score(y_true, y_pred)
        aris[method] = ari

    # Case 2: labels don't match → try all permutations
    else:
        unique_pred = y_pred.unique()
        unique_true = y_true.unique()

        best_ari = -1

        # Only permute if number of clusters matches
        if len(unique_pred) == len(unique_true):

            for perm in itertools.permutations(unique_true):
                mapping = dict(zip(unique_pred, perm))
                mapped_pred = y_pred.map(mapping)

                ari = adjusted_rand_score(y_true, mapped_pred)

                if ari > best_ari:
                    best_ari = ari

            aris[method] = best_ari

        else:
            # Direct ARI (label names don't matter mathematically)
            aris[method] = adjusted_rand_score(y_true, y_pred)

# -------------------------------------------------
# Print results
# -------------------------------------------------
print("Adjusted Rand Index per method:")
for k, v in aris.items():
    print(f"{k}: {v:.4f}")