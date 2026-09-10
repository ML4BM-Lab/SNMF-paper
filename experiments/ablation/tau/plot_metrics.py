
import sys
import os

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr
from skimage.metrics import structural_similarity as ssim

from scipy.stats import wilcoxon

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


def add_reference_significance(
    ax,
    data,
    x_col,
    y_col,
    order,
    reference_value,
    alternative,
):
    ref_vals = data[np.isclose(data[x_col], reference_value)][y_col].dropna()
    if ref_vals.empty:
        return

    ymin, ymax = ax.get_ylim()
    y_range = ymax - ymin
    if not np.isfinite(y_range) or y_range <= 0:
        y_range = max(abs(data[y_col].max()), 1.0)

    gap = y_range * 0.045
    max_annotation_y = ymax

    for value in order:
        if np.isclose(value, reference_value):
            continue

        vals = data[np.isclose(data[x_col], value)][y_col].dropna()
        if vals.empty:
            continue

        _, p = wilcoxon(vals, ref_vals, alternative=alternative)
        stars = pvalue_to_stars(p)
        y = vals.max() + gap

        ax.text(
            order.index(value),
            y,
            stars,
            ha="center",
            va="bottom",
            fontsize=18,
            zorder=10,
        )
        max_annotation_y = max(max_annotation_y, y)

    ax.set_ylim(ymin, max_annotation_y + gap)


sns.set_theme(style="whitegrid")
sns.set_context(
    "talk",
    font_scale=1.4,
    rc={
        "axes.titlesize": 25,
        "axes.labelsize": 22,
        "xtick.labelsize": 17,
        "ytick.labelsize": 17,
        "legend.fontsize": 17,
        "axes.linewidth": 1.4,
        "grid.linewidth": 0.7,
    },
)

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

order = sorted(rmse_df["value"].unique())
reference_value = 1.0
if not any(np.isclose(value, reference_value) for value in order):
    raise ValueError(
        f"Reference tau '{reference_value}' was not found in results: "
        + ", ".join(map(str, order))
    )

labels = [
    f"{value:g}" if value < 1 else f"{value:g} (NMF)"
    for value in order
]

palettes = {
    "rmse": sns.color_palette("Blues",  n_colors=len(order) + 2)[1:-1],
    "jsd":  sns.color_palette("Greens", n_colors=len(order) + 2)[1:-1],
    "pcc":  sns.color_palette("Purples", n_colors=len(order) + 2)[1:-1],
    "ssim": sns.color_palette("Oranges", n_colors=len(order) + 2)[1:-1],
}

# Plot boxplot
## RMSE
plt.figure(figsize=(11, 6))
ax = sns.boxplot(
    data=rmse_df,
    x="value",
    y="rmse",
    hue="value",
    order=order,
    hue_order=order,
    palette=palettes["rmse"],
    legend=False
)

add_reference_significance(
    ax,
    rmse_df,
    "value",
    "rmse",
    order,
    reference_value,
    alternative="less"
)

ax.set_xlabel("")
ax.set_ylabel("RMSE")

ax.set_xticks(range(len(order)))
ax.set_xticklabels(labels)

plt.tight_layout()
os.makedirs(os.path.join(results_path, "plots"), exist_ok=True)
plt.savefig(os.path.join(results_path, "plots", "rmse_comparison.png"), dpi=300)
plt.savefig(os.path.join(results_path, "plots", "rmse_comparison.pdf"), dpi=300)
plt.close()

## Jensen-Shannon divergence
plt.figure(figsize=(11, 6))
ax = sns.boxplot(
    data=jsd_df,
    x="value",
    y="jsd",
    hue="value",
    order=order,
    hue_order=order,
    palette=palettes["jsd"],
    legend=False
)

add_reference_significance(
    ax,
    jsd_df,
    "value",
    "jsd",
    order,
    reference_value,
    alternative="less"
)

ax.set_xlabel("")
ax.set_ylabel("JSD")

ax.set_xticks(range(len(order)))
ax.set_xticklabels(labels)

plt.tight_layout()
plt.savefig(os.path.join(results_path, "plots", "jsd_comparison.png"), dpi=300)
plt.savefig(os.path.join(results_path, "plots", "jsd_comparison.pdf"), dpi=300)
plt.close()

## PCC
plt.figure(figsize=(11, 6))
ax = sns.boxplot(
    data=pcc_df,
    x="value",
    y="pcc",
    hue="value",
    order=order,
    hue_order=order,
    palette=palettes["pcc"],
    legend=False
)

add_reference_significance(
    ax,
    pcc_df,
    "value",
    "pcc",
    order,
    reference_value,
    alternative="greater"
)

ax.set_xlabel("")
ax.set_ylabel("PCC")

ax.set_xticks(range(len(order)))
ax.set_xticklabels(labels)

plt.tight_layout()
plt.savefig(os.path.join(results_path, "plots", "pcc_comparison.png"), dpi=300)
plt.savefig(os.path.join(results_path, "plots", "pcc_comparison.pdf"), dpi=300)
plt.close()

## SSIM
plt.figure(figsize=(11, 6))
ax = sns.boxplot(
    data=ssim_df,
    x="value",
    y="ssim",
    hue="value",
    order=order,
    hue_order=order,
    palette=palettes["ssim"],
    legend=False
)

add_reference_significance(
    ax,
    ssim_df,
    "value",
    "ssim",
    order,
    reference_value,
    alternative="greater"
)

ax.set_xlabel("")
ax.set_ylabel("SSIM")

ax.set_xticks(range(len(order)))
ax.set_xticklabels(labels)

plt.tight_layout()
plt.savefig(os.path.join(results_path, "plots", "ssim_comparison.png"), dpi=300)
plt.savefig(os.path.join(results_path, "plots", "ssim_comparison.pdf"), dpi=300)
plt.close()
