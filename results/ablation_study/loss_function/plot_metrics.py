import sys
import os

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.spatial.distance import jensenshannon

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
    return "ns"


def add_significance(ax, data, x_col, y_col, order, pairs):
    ymax = data[y_col].max()
    y_offset = ymax * 0.05
    height = ymax * 0.02

    current_y = ymax + y_offset

    for g1, g2 in pairs:
        d1 = data[data[x_col] == g1][y_col]
        d2 = data[data[x_col] == g2][y_col]

        stat, p = mannwhitneyu(d1, d2, alternative="greater")
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

for dirpath, subdirs, files in os.walk(results_path):
    for f in files:
        if f == "SNMF_proportions.csv":
            value = dirpath.split("/")[-1]
            proportions = pd.read_csv(os.path.join(dirpath, f), index_col=0)

            proportions = proportions.loc[ground_truth.index]

            rmse = np.sqrt(((proportions - ground_truth) ** 2).mean(axis=1))
            rmses[value] = rmse.values

            P = proportions.values
            Q = ground_truth.values

            jsd = np.array([
                jensenshannon(p, q)**2
                for p, q in zip(P, Q)
            ])

            jsds[value] = jsd

# Convert to DataFrame (long format)
rmse_df = pd.DataFrame(
    [(val, v) for val, values in rmses.items() for v in values],
    columns=["value", "rmse"]
)

jsd_df = pd.DataFrame(
    [(val, v) for val, values in jsds.items() for v in values],
    columns=["value", "jsd"]
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