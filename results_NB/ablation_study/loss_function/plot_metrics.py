import sys
import os

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr
from skimage.metrics import structural_similarity as ssim

from scipy.stats import mannwhitneyu

def pvalue_to_stars(p):
    if p < 1e-4:
        return "****"
    elif p < 1e-3:
        return "***"
    elif p < 1e-2:
        return "**"
    elif p < 0.05:
        return "*"
    return "-"


def add_significance(ax, data, x_col, y_col, order, pairs, alternative="greater"):
    ymax = data[y_col].max()
    y_offset = ymax * 0.05
    height = ymax * 0.02

    current_y = ymax + y_offset

    for g1, g2 in pairs:
        d1 = data[data[x_col] == g1][y_col]
        d2 = data[data[x_col] == g2][y_col]

        _, p = mannwhitneyu(d1, d2, alternative=alternative)
        stars = pvalue_to_stars(p)

        x1 = order.index(g1)
        x2 = order.index(g2)

        ax.plot([x1, x1, x2, x2],
                [current_y, current_y+height, current_y+height, current_y],
                lw=1.5, c="black")

        ax.text((x1+x2)/2,
                current_y+height,
                stars,
                ha="center",
                va="bottom",
                fontsize=12)

        current_y += y_offset

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
            value = dirpath.split("/")[-1]
            proportions = pd.read_csv(os.path.join(dirpath, f), index_col=0)

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
                if np.std(p) != 0 and np.std(q) != 0
                else np.nan
                for p, q in zip(P, Q)
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

# ---------- FIXED ORDER ----------
order = sorted(rmse_df["value"].unique(), reverse=True)

# pairwise comparisons for significance
pairs = []
for i in range(len(order)):
    for j in range(i+1, len(order)):
        pairs.append((order[i], order[j]))

# ---------- RMSE PLOT ----------
plt.figure(figsize=(10,6))
ax = sns.violinplot(
    data=rmse_df,
    x="value",
    y="rmse",
    hue="value",
    palette=sns.color_palette("Set2"),
    order=order
)

add_significance(
    ax,
    rmse_df,
    "value",
    "rmse",
    order,
    pairs
)

ax.set_xlabel("Reconstruction")
ax.set_ylabel("RMSE")

plt.tight_layout()

os.makedirs(os.path.join(results_path, "plots"), exist_ok=True)

plt.savefig(os.path.join(results_path, "plots", "rmse_comparison.png"), dpi=300)
plt.savefig(os.path.join(results_path, "plots", "rmse_comparison.pdf"), dpi=300)
plt.close()


# ---------- JSD PLOT ----------
plt.figure(figsize=(10,6))
ax = sns.violinplot(
    data=jsd_df,
    x="value",
    y="jsd",
    hue="value",
    palette=sns.color_palette("Set2"),
    order=order
)

add_significance(
    ax,
    jsd_df,
    "value",
    "jsd",
    order,
    pairs
)

ax.set_xlabel("Reconstruction")
ax.set_ylabel("Jensen-Shannon divergence")

plt.tight_layout()

plt.savefig(os.path.join(results_path, "plots", "jsd_comparison.png"), dpi=300)
plt.savefig(os.path.join(results_path, "plots", "jsd_comparison.pdf"), dpi=300)
plt.close()

# ---------- PCC PLOT ----------
plt.figure(figsize=(10,6))
ax = sns.violinplot(
    data=pcc_df,
    x="value",
    y="pcc",
    hue="value",
    palette=sns.color_palette("Set2"),
    order=order
)

add_significance(
    ax,
    pcc_df,
    "value",
    "pcc",
    order,
    pairs
)

ax.set_xlabel("Reconstruction")
ax.set_ylabel("Pearson Correlation Coefficient")

plt.tight_layout()

plt.savefig(os.path.join(results_path, "plots", "pcc_comparison.png"), dpi=300)
plt.savefig(os.path.join(results_path, "plots", "pcc_comparison.pdf"), dpi=300)
plt.close()

# ---------- SSIM PLOT ----------
plt.figure(figsize=(10,6))
ax = sns.violinplot(
    data=ssim_df,
    x="value",
    y="ssim",
    hue="value",
    palette=sns.color_palette("Set2"),
    order=order
)

add_significance(
    ax,
    ssim_df,
    "value",
    "ssim",
    order,
    pairs
)

ax.set_xlabel("Reconstruction")
ax.set_ylabel("SSIM")

plt.tight_layout()

plt.savefig(os.path.join(results_path, "plots", "ssim_comparison.png"), dpi=300)
plt.savefig(os.path.join(results_path, "plots", "ssim_comparison.pdf"), dpi=300)
plt.close()