import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
from scipy.stats import ttest_ind

# === Configuration ===
MEM_CONVERSION = {'T': 2, 'G': 1, 'M': 0, 'K': -1}
METHOD_ORDER = ["SMART", "STdeconvolve", "starfysh", "CARD", "RETROFIT", "BayesTME", "SpiceMix", "SNMF", "SNMF_v2"]

# Global seaborn style
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


# === Helper functions ===
def get_csv_files(folder, hungarian):
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


def rmse(output, ground_truth):
    return np.sqrt(((output - ground_truth) ** 2).mean(axis=1)).tolist()


def get_metrics(output, fpath, ground_truth):
    mem, time = parse_sacct_output(fpath)
    return {'time': [time], 'mem': [mem], 'rmse': rmse(output, ground_truth)}


def plot_metric_bar(df, y_label, title, filename, log_scale=False):
    plt.figure(figsize=(9, 6))
    methods_present = [m for m in METHOD_ORDER if m in df["Method"].unique()]
    palette = [METHOD_COLOR_MAP[m] for m in methods_present]

    ax = sns.barplot(
        data=df, x="CellTypes", y=y_label,
        hue="Method", hue_order=methods_present,
        palette=palette, edgecolor='black', linewidth=1.3
    )

    if log_scale:
        ax.set_yscale("log")

    # ax.set_title(title, pad=15, weight='bold', fontsize=18)
    ax.set_ylabel(y_label)
    ax.set_xlabel("Number of Cell Types")

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
    plt.savefig(filename, dpi=300)
    plt.close()


def aggregate_metrics(experiments):
    """Combine metrics from multiple experiments into one DataFrame per metric."""
    rmse_data, time_data, mem_data = [], [], []

    for n_cells, metrics in experiments.items():
        for m, vals in metrics.items():
            # RMSE (multiple per run)
            for val in vals["rmse"]:
                rmse_data.append({"CellTypes": int(n_cells), "Method": m, "RMSE": val})
            # Single values for time/mem
            time_data.append({"CellTypes": int(n_cells), "Method": m, "Time (s)": np.nanmean(vals["time"])})
            mem_data.append({"CellTypes": int(n_cells), "Method": m, "Memory (MB)": np.nanmean(vals["mem"])})

    return (
        pd.DataFrame(rmse_data),
        pd.DataFrame(time_data),
        pd.DataFrame(mem_data)
    )


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <root_folder> <hungarian:true/false>")
        sys.exit(1)

    root_folder = sys.argv[1]
    ground_truth_folder = sys.argv[2]
    hungarian = sys.argv[3].lower() == "true"

    experiments = {}

    for n_cells in ["5", "8", "10", "15"]:
        exp_folder = os.path.join(root_folder, n_cells)
        gt_path = os.path.join(ground_truth_folder, f"k{n_cells}", f"proportions_k{n_cells}.csv")

        if not os.path.isdir(exp_folder):
            print(f"⚠️ Skipping missing folder: {exp_folder}")
            continue
        if not os.path.exists(gt_path):
            print(f"⚠️ Missing ground truth for {n_cells} cell types: {gt_path}")
            continue

        print(f"📂 Processing {n_cells} cell types...")
        ground_truth = pd.read_csv(gt_path, index_col=0)
        csv_files = get_csv_files(exp_folder, hungarian)

        metrics = {}
        for fpath in csv_files:
            parent_folder = os.path.basename(os.path.dirname(fpath))
            output = pd.read_csv(fpath, index_col=0)
            if output.shape == ground_truth.shape:
                metrics[parent_folder] = get_metrics(output, fpath, ground_truth)
        experiments[n_cells] = metrics

    # Aggregate across experiments
    rmse_df, time_df, mem_df = aggregate_metrics(experiments)

    plots_dir = os.path.join(root_folder, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    plot_metric_bar(rmse_df, "RMSE", "RMSE vs. Number of Cell Types", os.path.join(plots_dir, "rmse_vs_celltypes.png"))
    plot_metric_bar(time_df, "Time (s)", "Execution Time vs. Number of Cell Types",
                    os.path.join(plots_dir, "time_vs_celltypes.png"), log_scale=True)
    plot_metric_bar(mem_df, "Memory (MB)", "Memory Usage vs. Number of Cell Types",
                    os.path.join(plots_dir, "memory_vs_celltypes.png"))

    print("✅ Plots saved in 'plots/' folder under root directory.")


if __name__ == "__main__":
    main()
