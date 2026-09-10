import sys
import os

import scanpy as sc
import pandas as pd

from sklearn.metrics import confusion_matrix
from scipy.optimize import linear_sum_assignment
import numpy as np

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
ax = sc.pl.spatial(
    adata,
    color="Region",
    show=False,
    frameon=False,
    na_in_legend=False,
    palette="tab20",
    title="",
)
ax[0].set_xlabel("")
ax[0].set_ylabel("")

plt.savefig(os.path.join(plots_dir, "ground_truth.png"),
            dpi=300,
            bbox_inches="tight")
plt.close()

# Store ground-truth color mapping
gt_categories = adata.obs["Region"].cat.categories
gt_colors = dict(zip(gt_categories, adata.uns["Region_colors"]))

for method in methods:

    pred_key = f"{method}_major_ct"
    
    # Remove NA
    mask = adata.obs[pred_key].notna() & adata.obs["Region"].notna()
    y_true = adata.obs.loc[mask, "Region"]
    y_pred = adata.obs.loc[mask, pred_key]

    # Build confusion matrix
    true_labels = np.unique(y_true)
    pred_labels = np.unique(y_pred)

    true_label_map = {lab: i for i, lab in enumerate(true_labels)}
    y_true_cat = np.array([true_label_map[x] for x in y_true])

    pred_label_map = {lab: i for i, lab in enumerate(pred_labels)}
    y_pred_cat = np.array([pred_label_map[x] for x in y_pred])

    cm = confusion_matrix(y_true_cat, y_pred_cat)

    # Hungarian matching (maximize overlap)
    row_ind, col_ind = linear_sum_assignment(-cm)

    # Create mapping
    mapping = {
        pred_labels[j]: true_labels[i]
        for i, j in zip(row_ind, col_ind)
        if j < len(pred_labels)
    }

    # Apply remapping
    remapped = adata.obs[pred_key].map(mapping)
    adata.obs[f"{method}_aligned"] = remapped.astype("category")

    # Force same category order as ground truth
    adata.obs[f"{method}_aligned"] = (
        adata.obs[f"{method}_aligned"]
        .cat.set_categories(gt_categories)
    )

    # Assign same colors
    adata.uns[f"{method}_aligned_colors"] = [
        gt_colors[c] for c in gt_categories
    ]

    # Plot
    ax = sc.pl.spatial(
        adata,
        color=f"{method}_aligned",
        show=False,
        frameon=False,
        na_in_legend=False,
        legend_loc="none",
        title=""
    )

    ax[0].set_xlabel("")
    ax[0].set_ylabel("")

    plt.savefig(
        os.path.join(plots_dir, f"{method}.png"),
        dpi=300,
        bbox_inches="tight"
    )
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

    # Direct ARI (label names don't matter mathematically)
    aris[method] = adjusted_rand_score(y_true, y_pred)

# -------------------------------------------------
# Print and save results
# -------------------------------------------------
with open(os.path.join(results_path, "ari.txt"), 'w') as f:
    f.write("Adjusted Rand Index per method:\n")
    for k, v in aris.items():
        line = f"{k}: {v:.4f}\n"
        print(line.strip())
        f.write(line)