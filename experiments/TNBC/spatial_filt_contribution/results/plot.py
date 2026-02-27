import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================
# Setup
# ==========================
results_path = sys.argv[1]
proportions_path = sys.argv[2]

ground_truth = pd.read_csv(proportions_path, index_col=0)
results = {}

def rmse(output):
    return np.sqrt(((output - ground_truth) ** 2).mean(axis=1)).tolist()

# ==========================
# Collect results
# ==========================
for root, _, files in os.walk(results_path):
    for f in files:
        if f in ("SNMF_proportions.csv", "hungarian_proportions.csv"):
            method = root.split("/")[-1]
            n_runs = root.split("/")[-2]
            output = pd.read_csv(os.path.join(root, f), index_col=0)
            results.setdefault(n_runs, {})[method] = rmse(output)

# ==========================
# Prepare long-format DataFrame
# ==========================
plot_data = []
for n_run, methods in results.items():
    for method, values in methods.items():
        for v in values:
            plot_data.append({
                "n_run": int(n_run.split("_")[-1]),
                "method": method.replace("_proportions", ""),
                "rmse": v
            })

df_plot = pd.DataFrame(plot_data).sort_values("n_run")

# ==========================
# Style configuration
# ==========================
sns.set_theme(style="whitegrid", context="paper", font_scale=0.9)
method_order = ["NMF", "SNMF"]
set2_colors = sns.color_palette("Set2", n_colors=2)
palette = dict(zip(method_order, set2_colors))

# ==========================
# --- Boxplot ---
# ==========================
plt.figure(figsize=(7, 4))
sns.boxplot(
    data=df_plot,
    x="n_run",
    y="rmse",
    hue="method",
    hue_order=method_order,
    palette=palette,
    width=0.6,
    fliersize=2,
    linewidth=1.0
)

plt.xlabel("Number of Runs", labelpad=6)
plt.ylabel("RMSE", labelpad=6)
plt.title("")  # no title
sns.despine(offset=5, trim=True)

# legend closer to plot, no title
plt.legend(
    title=None,
    frameon=False,
    bbox_to_anchor=(1.01, 0.5),
    loc="center left"
)

plt.tight_layout(rect=[0, 0, 0.93, 1])  # less right margin
plt.savefig(os.path.join(results_path, "box_plots.png"), dpi=300, transparent=False)
plt.close()

# ==========================
# --- Barplot ---
# ==========================
plot_data = []
for n_run, methods in results.items():
    for method, values in methods.items():
        for v in values:  # keep each RMSE value
            plot_data.append({
                "n_run": int(n_run.split("_")[-1]),
                "method": method.replace("_proportions", ""),
                "rmse": v
            })

df_plot = pd.DataFrame(plot_data).sort_values("n_run")

plt.figure(figsize=(7, 4))
ax = sns.barplot(
    data=df_plot,
    x="n_run",
    y="rmse",
    hue="method",
    hue_order=method_order,
    palette=palette,
    capsize=0.15,
    errcolor="gray",
    errwidth=1.0
)

plt.xlabel("Number of Runs", labelpad=6)
plt.ylabel("Mean RMSE", labelpad=6)
plt.title("")  # no title
sns.despine(offset=5, trim=True)

plt.legend(
    title=None,
    frameon=False,
    bbox_to_anchor=(1.01, 0.5),
    loc="center left"
)

plt.tight_layout(rect=[0, 0, 0.93, 1])
plt.savefig(os.path.join(results_path, "bar_plots.png"), dpi=300, transparent=False)
plt.close()
