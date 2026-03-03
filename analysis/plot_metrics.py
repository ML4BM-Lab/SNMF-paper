import os
import sys

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns

from scipy.spatial.distance import jensenshannon

from scipy.stats import ttest_ind

# === Configuration ===
MEM_CONVERSION = {'T': 2, 'G': 1, 'M': 0, 'K': -1}
METHOD_ORDER = [
    "SNMF",
    "STdeconvolve",
    "SpiceMix",
    "CARD",
    "SMART",
    "RETROFIT",
    "starfysh",
    "BayesTME"
]

# Global seaborn style with fine-tuned fonts
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

METHOD_COLORS = sns.color_palette("Spectral", n_colors=len(METHOD_ORDER))
METHOD_COLOR_MAP = dict(zip(METHOD_ORDER, METHOD_COLORS))
# METHOD_COLOR_MAP["CARD"] = "#FFC84D"

ground_truth = None
folder = None
hungarian = None


# === Helper functions ===
def get_csv_files():
    global folder, hungarian
    csv_files = []
    for root, _, files in os.walk(folder):
        for f in files:
            if f.lower().endswith(".csv") and "gene_sig" not in f and (('hungarian' in f) == hungarian):
                csv_files.append(os.path.join(root, f))
    return csv_files


def parse_sacct_output(fpath):
    mem, time = None, None
    sacct_file = os.path.join(os.path.dirname(fpath), "sacct.log")
    if not os.path.exists(sacct_file):
        return np.nan, np.nan
    with open(sacct_file, 'r') as f:
        for line in f.readlines():
            if "batch" in line and "COMPLETED" in line:
                mem, time = line.split()[2:4]
    if mem is not None:
        mem = float(mem[:-1]) * 1024 ** MEM_CONVERSION[mem[-1]]  # Convert to MB
    if time is not None:
        h, m, s = map(int, time.split(":"))
        time = h * 3600 + m * 60 + s
    return mem, time


def rmse(output):
    global ground_truth
    return np.sqrt(((output - ground_truth) ** 2).mean(axis=1)).tolist()

def jsd(output):
    global ground_truth
    P = output.values
    Q = ground_truth.values
    return [jensenshannon(p, q)**2 for p, q in zip(P, Q)]

def get_metrics(output, fpath):
    mem, time = parse_sacct_output(fpath)
    return {
        'time': [time],
        'mem': [mem],
        'rmse': rmse(output),
        'jsd': jsd(output)
    }


def significance_stars(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "-"


def add_rmse_significance_stars(ax, data, methods_present):
    """Add significance stars above bars."""
    if "SNMF" not in data:
        return

    ref_vals = np.array(data["SNMF"]["rmse"])
    y_max_all = 0

    for i, method in enumerate(methods_present):
        if method == "SNMF":
            continue
        vals = np.array(data[method]["rmse"])
        _, pval = ttest_ind(vals, ref_vals, equal_var=False, nan_policy='omit')
        star = significance_stars(pval)

        y_max = np.nanmax(vals)
        y = y_max * 1.05
        y_max_all = max(y_max_all, y)
        ax.text(i, y, star, ha='center', va='bottom',
                fontsize=14, fontweight='bold', color='black')

    ylim = ax.get_ylim()
    ax.set_ylim(ylim[0], max(ylim[1], y_max_all * 1.15))


def plot_bar(data, y_label, title, filename, log_scale=False):
    global folder
    plt.figure(figsize=(9.5, 5))
    methods_present = [m for m in METHOD_ORDER if m in data["Method"].unique()]
    palette = [METHOD_COLOR_MAP[m] for m in methods_present]

    ax = sns.barplot(
        data=data, x="Method", y=y_label,
        order=methods_present, palette=palette,
        edgecolor='black', linewidth=1.3
    )

    if log_scale:
        ax.set_yscale("log")

    ax.set_xlabel("")
    ax.set_ylabel(y_label, fontsize=15)
    ax.set_title(title, pad=15, weight='bold', fontsize=18)

    ax.set_xticks([])

    legend_handles = [
        Patch(facecolor=METHOD_COLOR_MAP[m], edgecolor='black', label=m)
        for m in methods_present
    ]
    ax.legend(
        handles=legend_handles,
        title="Method", title_fontsize=13,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=True
    )

    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(folder, "plots", filename), dpi=300)
    plt.savefig(os.path.join(folder, "plots", f"{filename.split('.')[0]}.pdf"), dpi=300)
    plt.close()


def plot_rmse(metrics):
    global folder
    methods_present = [m for m in METHOD_ORDER if m in metrics]
    palette = [METHOD_COLOR_MAP[m] for m in methods_present]

    rmse_data = []
    for m in methods_present:
        for val in metrics[m]['rmse']:
            rmse_data.append({'Method': m, 'RMSE': val})
    rmse_df = pd.DataFrame(rmse_data)

    plt.figure(figsize=(11, 6))
    ax = sns.boxplot(
        x="Method", y="RMSE", data=rmse_df,
        order=methods_present, palette=palette,
        width=0.65, fliersize=2, linewidth=1.2
    )

    add_rmse_significance_stars(ax, metrics, methods_present)
    ax.set_title("RMSE Distribution per Method", pad=15, weight='bold', fontsize=18)
    ax.set_xlabel("")
    ax.set_ylabel("RMSE", fontsize=15)
    ax.set_xticks([])

    legend_handles = [
        Patch(facecolor=METHOD_COLOR_MAP[m], edgecolor='black', label=m)
        for m in methods_present
    ]
    ax.legend(
        handles=legend_handles,
        title="Method", title_fontsize=13,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=True
    )

    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(folder, "plots", "rmse_comparison.png"), dpi=300)
    plt.savefig(os.path.join(folder, "plots", "rmse_comparison.pdf"), dpi=300)
    plt.close()

def plot_jsd(metrics):
    global folder

    methods_present = [m for m in METHOD_ORDER if m in metrics]
    palette = [METHOD_COLOR_MAP[m] for m in methods_present]

    jsd_data = []
    for m in methods_present:
        for val in metrics[m]['jsd']:
            jsd_data.append({'Method': m, 'JSD': val})

    jsd_df = pd.DataFrame(jsd_data)

    plt.figure(figsize=(11, 6))
    ax = sns.boxplot(
        x="Method", y="JSD", data=jsd_df,
        order=methods_present, palette=palette,
        width=0.65, fliersize=2, linewidth=1.2
    )

    # Significance vs SNMF
    if "SNMF" in metrics:
        ref_vals = np.array(metrics["SNMF"]["jsd"])
        y_max_all = 0

        for i, method in enumerate(methods_present):
            if method == "SNMF":
                continue
            vals = np.array(metrics[method]["jsd"])
            _, pval = ttest_ind(vals, ref_vals, equal_var=False, nan_policy='omit')
            star = significance_stars(pval)

            y_max = np.nanmax(vals)
            y = y_max * 1.05
            y_max_all = max(y_max_all, y)
            ax.text(i, y, star, ha='center', va='bottom',
                    fontsize=14, fontweight='bold', color='black')

        ylim = ax.get_ylim()
        ax.set_ylim(ylim[0], max(ylim[1], y_max_all * 1.15))

    ax.set_title("Jensen-Shannon Divergence per Method", pad=15, weight='bold', fontsize=18)
    ax.set_xlabel("")
    ax.set_ylabel("Jensen-Shannon Divergence", fontsize=15)
    ax.set_xticks([])

    legend_handles = [
        Patch(facecolor=METHOD_COLOR_MAP[m], edgecolor='black', label=m)
        for m in methods_present
    ]

    ax.legend(
        handles=legend_handles,
        title="Method", title_fontsize=13,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=True
    )

    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(folder, "plots", "jsd_comparison.png"), dpi=300)
    plt.savefig(os.path.join(folder, "plots", "jsd_comparison.pdf"), dpi=300)
    plt.close()


def plot_metrics(metrics):
    global folder
    os.makedirs(os.path.join(folder, "plots"), exist_ok=True)
    methods_present = [m for m in METHOD_ORDER if m in metrics]

    time_data = pd.DataFrame([
        {"Method": m, "Time (s)": np.nanmean(metrics[m]['time'])} for m in methods_present
    ])
    plot_bar(time_data, "Time (s)", "Execution Time Comparison", "time_comparison.png", log_scale=True)

    mem_data = pd.DataFrame([
        {"Method": m, "Memory (MB)": np.nanmean(metrics[m]['mem'])} for m in methods_present
    ])
    plot_bar(mem_data, "Memory (MB)", "Memory Usage Comparison", "memory_comparison.png", log_scale=False)

    plot_rmse(metrics)
    plot_jsd(metrics)


def main():
    global ground_truth, folder, hungarian

    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} <folder_path> <ground_truth_path> <hungarian:true/false>")
        sys.exit(1)

    folder = sys.argv[1]
    ground_truth = pd.read_csv(sys.argv[2], index_col=0)
    hungarian = sys.argv[3].lower() == "true"

    csv_files = get_csv_files()
    if not csv_files:
        print("No CSV files found.")
        sys.exit(0)

    metrics = {}
    for fpath in csv_files:
        parent_folder = os.path.basename(os.path.dirname(fpath))
        output = pd.read_csv(fpath, index_col=0)
        if output.shape == ground_truth.shape:
            metrics[parent_folder] = get_metrics(output, fpath)

    plot_metrics(metrics)
    print("✅ Plots saved in 'plots/' folder.")


if __name__ == "__main__":
    main()
