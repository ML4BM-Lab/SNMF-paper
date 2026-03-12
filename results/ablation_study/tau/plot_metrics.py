import sys
import os

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.spatial.distance import jensenshannon

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
            value = float(dirpath.split("/")[-1][1:])
            proportions = pd.read_csv(os.path.join(dirpath, f), index_col=0)

            # Ensure same order
            proportions = proportions.loc[ground_truth.index]

            rmse = np.sqrt(((proportions - ground_truth) ** 2).mean(axis=1))
            rmses[value] = rmse.values

            # Jensen-Shannon divergence
            P = proportions.values
            Q = ground_truth.values

            # Compute JSD per row (scipy returns sqrt(JS), so square it)
            jsd = np.array([
                jensenshannon(p, q)**2
                for p, q in zip(P, Q)
            ])

            jsds[value] = jsd

# Convert to DataFrame (long format)
rmse_df = pd.DataFrame(
    [(val, value) for val, rmses in rmses.items() for value in rmses],
    columns=["value", "rmse"]
)

jsd_df = pd.DataFrame(
    [(val, value) for val, jsd in jsds.items() for value in jsd],
    columns=["value", "jsd"]
)

# Sort by value
rmse_df = rmse_df.sort_values("value")
jsd_df = jsd_df.sort_values("value")

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