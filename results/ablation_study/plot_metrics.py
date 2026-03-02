import sys
import os

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

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

results = {}

for dirpath, subdirs, files in os.walk(results_path):
    for f in files:
        if f == "SNMF_proportions.csv":
            value = float(dirpath.split("/")[-1][1:])
            proportions = pd.read_csv(os.path.join(dirpath, f), index_col=0)

            # Ensure same order
            proportions = proportions.loc[ground_truth.index]

            rmse = np.sqrt(((proportions - ground_truth) ** 2).mean(axis=1))
            results[value] = rmse.values

# Convert to DataFrame (long format)
rmse_df = pd.DataFrame(
    [(val, value) for val, rmses in results.items() for value in rmses],
    columns=["value", "rmse"]
)

# Sort by value
rmse_df = rmse_df.sort_values("value")

# Plot boxplot
plt.figure(figsize=(10, 6))
ax = sns.boxplot(data=rmse_df, x="value", y="rmse", hue="value", legend=False)

ax.set_xlabel("Tau")
ax.set_ylabel("RMSE")

labels = [f"{g}" if g < 7 else f"{g} (NMF)" for g in rmse_df["value"].unique()]
ax.set_xticklabels(labels)

plt.tight_layout()
os.makedirs(os.path.join(results_path, "plots"), exist_ok=True)
plt.savefig(os.path.join(results_path, "plots", "rmse_comparison.png"), dpi=300)
plt.show()