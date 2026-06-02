
import sys
import os

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr
from skimage.metrics import structural_similarity as ssim

sns.set_theme(style="whitegrid")
sns.set_context("talk", font_scale=1.1, rc={
    "axes.titlesize": 18,
    "axes.labelsize": 15,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "axes.linewidth": 1.2,
    "grid.linewidth": 0.6
})

results_path = sys.argv[1]
ground_truth = pd.read_csv(sys.argv[2], index_col=0)

rmses = {}
jsds = {}
pccs = {}
ssims = {}

for dirpath, subdirs, files in os.walk(results_path):
    for f in files:
        if f == "SNMF_proportions.csv":
            value = float(dirpath.split("/")[-1][1:])
            proportions = pd.read_csv(os.path.join(dirpath, f), index_col=0)

            # Ensure same order
            proportions = proportions.loc[ground_truth.index]

            # RMSE
            rmse = np.sqrt(((proportions - ground_truth) ** 2).mean(axis=1))
            rmses[value] = rmse.values

            # JSD
            P = proportions.values
            Q = ground_truth.values

            jsd = np.array([
                jensenshannon(p, q)**2
                for p, q in zip(P, Q)
            ])

            jsds[value] = jsd

            # PCC
            pcc = np.array([
                pearsonr(p, q)[0]
                for p, q in zip(P, Q)
                if np.std(p) != 0 and np.std(q) != 0
            ])
            pccs[value] = pcc

            # SSIM
            coords = np.array([
                list(map(int, idx.split("x")))
                for idx in proportions.index
            ])

            rows = coords[:, 0]
            cols = coords[:, 1]

            H = rows.max() + 1
            W = cols.max() + 1

            vals = []

            for celltype in ground_truth.columns:

                pred_img = np.zeros((H, W), dtype=float)
                true_img = np.zeros((H, W), dtype=float)

                pred_vals = proportions[celltype.replace("-", ".").replace("&", ".")].values
                true_vals = ground_truth.loc[[f"{r}x{c}" for r,c in zip(rows, cols)], celltype].values

                pred_img[rows, cols] = pred_vals
                true_img[rows, cols] = true_vals

                data_range = max(
                    pred_img.max() - pred_img.min(),
                    true_img.max() - true_img.min(),
                    1e-8
                )

                score = ssim(
                    pred_img,
                    true_img,
                    data_range=data_range
                )

                vals.append(score)

            ssims[value] = np.array(vals)

# Convert to DataFrame (long format)
rmse_df = pd.DataFrame(
    [(val, v) for val, values in rmses.items() for v in values],
    columns=["value", "rmse"]
)

jsd_df = pd.DataFrame(
    [(val, v) for val, values in jsds.items() for v in values],
    columns=["value", "jsd"]
)

pcc_df = pd.DataFrame(
    [(val, v) for val, values in pccs.items() for v in values],
    columns=["value", "pcc"]
)

ssim_df = pd.DataFrame(
    [(val, v) for val, values in ssims.items() for v in values],
    columns=["value", "ssim"]
)

# Sort by value
rmse_df = rmse_df.sort_values("value")
jsd_df = jsd_df.sort_values("value")
pcc_df = pcc_df.sort_values("value")
ssim_df = ssim_df.sort_values("value")

# Plot boxplot
## RMSE
plt.figure(figsize=(10, 6))
ax = sns.boxplot(data=rmse_df, x="value", y="rmse", hue="value", legend=False)

ax.set_xlabel("Tau")
ax.set_ylabel("RMSE")

labels = [f"{g}" if g < 1 else f"{g} (NMF)" for g in rmse_df["value"].unique()]
ax.set_xticklabels(labels)

plt.tight_layout()
os.makedirs(os.path.join(results_path, "plots"), exist_ok=True)
plt.savefig(os.path.join(results_path, "plots", "rmse_comparison.png"), dpi=300)
plt.savefig(os.path.join(results_path, "plots", "rmse_comparison.pdf"), dpi=300)
plt.close()

## Jensen-Shannon divergence
plt.figure(figsize=(10, 6))
ax = sns.boxplot(data=jsd_df, x="value", y="jsd", hue="value", legend=False)

ax.set_xlabel("Tau")
ax.set_ylabel("Jensen-Shannon divergence")

labels = [f"{g}" if g < 1 else f"{g} (NMF)" for g in rmse_df["value"].unique()]
ax.set_xticklabels(labels)

plt.tight_layout()
plt.savefig(os.path.join(results_path, "plots", "jsd_comparison.png"), dpi=300)
plt.savefig(os.path.join(results_path, "plots", "jsd_comparison.pdf"), dpi=300)
plt.close()

## PCC
plt.figure(figsize=(10, 6))
ax = sns.boxplot(data=pcc_df, x="value", y="pcc", hue="value", legend=False)

ax.set_xlabel("Tau")
ax.set_ylabel("PCC")

labels = [f"{g}" if g < 1 else f"{g} (NMF)" for g in rmse_df["value"].unique()]
ax.set_xticklabels(labels)

plt.tight_layout()
plt.savefig(os.path.join(results_path, "plots", "pcc_comparison.png"), dpi=300)
plt.savefig(os.path.join(results_path, "plots", "pcc_comparison.pdf"), dpi=300)
plt.close()

## SSIM
plt.figure(figsize=(10, 6))
ax = sns.boxplot(data=ssim_df, x="value", y="ssim", hue="value", legend=False)

ax.set_xlabel("Tau")
ax.set_ylabel("SSIM")

labels = [f"{g}" if g < 1 else f"{g} (NMF)" for g in rmse_df["value"].unique()]
ax.set_xticklabels(labels)

plt.tight_layout()
plt.savefig(os.path.join(results_path, "plots", "ssim_comparison.png"), dpi=300)
plt.savefig(os.path.join(results_path, "plots", "ssim_comparison.pdf"), dpi=300)
plt.close()