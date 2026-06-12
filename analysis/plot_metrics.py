import os
import sys

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns

from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr
from skimage.metrics import structural_similarity as ssim

from scipy.stats import wilcoxon

# === Configuration ===
MEM_CONVERSION = {'T': 2, 'G': 1, 'M': 0, 'K': -1}
METHOD_ORDER = [
    "SNMF",
    "NMF",
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

def pcc(output):
    global ground_truth

    P = output.values
    Q = ground_truth.values

    vals = []

    for p, q in zip(P, Q):
        if np.std(p) != 0 and np.std(q) != 0:
            vals.append(pearsonr(p, q)[0])

    return vals

def compute_ssim(output):
    global ground_truth

    coords = np.array([
        list(map(int, idx.split("x")))
        for idx in output.index
    ])

    rows = coords[:, 0]
    cols = coords[:, 1]

    H = rows.max() + 1
    W = cols.max() + 1

    vals = []

    # ----------------------------------------
    # Compute SSIM per cell type
    # ----------------------------------------

    for celltype in ground_truth.columns:

        pred_img = np.zeros((H, W), dtype=float)
        true_img = np.zeros((H, W), dtype=float)

        pred_vals = output[celltype.replace("-", ".").replace("&", ".")].values
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

    return vals

def get_metrics(output, fpath):
    mem, time = parse_sacct_output(fpath)
    return {
        'time': [time],
        'mem': [mem],
        'rmse': rmse(output),
        'jsd': jsd(output),
        'pcc': pcc(output),
        'ssim': compute_ssim(output)
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


def add_significance_stars(ax, data, methods_present, metric="rmse", alternative='greater'):
    """Add significance stars above bars."""

    if "SNMF" not in data:
        return

    ref_vals = np.array(data["SNMF"][metric])
    y_max_all = 0

    for method in data.keys():
        if method == "SNMF":
            continue
        vals = np.array(data[method][metric])
        pval = wilcoxon(vals, ref_vals, alternative=alternative).pvalue # Wilcoxon signed-rank (not normal assumption)

        star = significance_stars(pval)

        y_max = np.nanmax(vals)
        y = y_max * 1.05
        y_max_all = max(y_max_all, y)
        ax.text(methods_present.index(method), y, star, ha='center', va='bottom',
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

def plot_distribution_metric(metrics, metric_name, ylabel, filename, alternative):
    global folder

    methods_present = [m for m in METHOD_ORDER if m in metrics]
    palette = [METHOD_COLOR_MAP[m] for m in methods_present]

    data_rows = []

    for m in methods_present:
        for val in metrics[m][metric_name]:
            data_rows.append({
                'Method': m,
                ylabel: val
            })

    df = pd.DataFrame(data_rows)

    plt.figure(figsize=(11, 6))

    ax = sns.boxplot(
        x='Method',
        y=ylabel,
        data=df,
        order=methods_present,
        palette=palette,
        width=0.65,
        fliersize=2,
        linewidth=1.2
    )

    add_significance_stars(ax, metrics, methods_present, metric_name, alternative=alternative)
    ax.set_title(f'{ylabel} per Method', pad=15, weight='bold', fontsize=18)
    ax.set_xlabel('')
    ax.set_ylabel(ylabel, fontsize=15)
    ax.set_xticks([])

    legend_handles = [
        Patch(facecolor=METHOD_COLOR_MAP[m], edgecolor='black', label=m)
        for m in methods_present
    ]

    ax.legend(
        handles=legend_handles,
        title='Method',
        title_fontsize=13,
        loc='center left',
        bbox_to_anchor=(1.02, 0.5),
        frameon=True
    )

    sns.despine()
    plt.tight_layout()

    plt.savefig(os.path.join(folder, 'plots', filename), dpi=300)
    plt.savefig(
        os.path.join(folder, 'plots', f"{filename.split('.')[0]}.pdf"),
        dpi=300
    )

    plt.close()


def plot_metrics(metrics):
    global folder
    os.makedirs(os.path.join(folder, "plots"), exist_ok=True)

    save_metric_summary(metrics)

    
    methods_present = [m for m in METHOD_ORDER if m in metrics]

    time_data = pd.DataFrame([
        {"Method": m, "Time (s)": np.nanmean(metrics[m]['time'])} for m in methods_present
    ])
    plot_bar(time_data, "Time (s)", "Execution Time Comparison", "time_comparison.png", log_scale=True)

    mem_data = pd.DataFrame([
        {"Method": m, "Memory (MB)": np.nanmean(metrics[m]['mem'])} for m in methods_present
    ])
    plot_bar(mem_data, "Memory (MB)", "Memory Usage Comparison", "memory_comparison.png", log_scale=False)

    plot_distribution_metric(
        metrics,
        metric_name='rmse',
        ylabel='RMSE',
        filename='rmse_comparison.png',
        alternative='greater',
    )
    plot_distribution_metric(
        metrics,
        metric_name='jsd',
        ylabel='JSD',
        filename='jsd_comparison.png',
        alternative='greater',
    )
    plot_distribution_metric(
        metrics,
        metric_name='pcc',
        ylabel='PCC',
        filename='pcc_comparison.png',
        alternative='less',
    )
    plot_distribution_metric(
        metrics,
        metric_name='ssim',
        ylabel='SSIM',
        filename='ssim_comparison.png',
        alternative='less',
    )


def save_metric_summary(metrics):
    global folder

    rows = []

    for method in [m for m in METHOD_ORDER if m in metrics]:
        rmse_vals = np.array(metrics[method]["rmse"])
        jsd_vals = np.array(metrics[method]["jsd"])
        pcc_vals = np.array(metrics[method]["pcc"])
        ssim_vals = np.array(metrics[method]["ssim"])

        rows.append({
            "Method": method,
            "Time": metrics[method]["time"][0],
            "Memory": metrics[method]["mem"][0],

            "RMSE_mean": np.nanmean(rmse_vals),
            "RMSE_var": np.nanvar(rmse_vals),
            "RMSE_Q1": np.nanpercentile(rmse_vals, 25),
            "RMSE_Q3": np.nanpercentile(rmse_vals, 75),
            "RMSE_IQR": np.nanpercentile(rmse_vals, 75) - np.nanpercentile(rmse_vals, 25),

            "JSD_mean": np.nanmean(jsd_vals),
            "JSD_var": np.nanvar(jsd_vals),
            "JSD_Q1": np.nanpercentile(jsd_vals, 25),
            "JSD_Q3": np.nanpercentile(jsd_vals, 75),
            "JSD_IQR": np.nanpercentile(jsd_vals, 75) - np.nanpercentile(jsd_vals, 25),

            "PCC_mean": np.nanmean(pcc_vals),
            "PCC_var": np.nanvar(pcc_vals),
            "PCC_Q1": np.nanpercentile(pcc_vals, 25),
            "PCC_Q3": np.nanpercentile(pcc_vals, 75),
            "PCC_IQR": np.nanpercentile(pcc_vals, 75) - np.nanpercentile(pcc_vals, 25),

            "SSIM_mean": np.nanmean(ssim_vals),
            "SSIM_var": np.nanvar(ssim_vals),
            "SSIM_Q1": np.nanpercentile(ssim_vals, 25),
            "SSIM_Q3": np.nanpercentile(ssim_vals, 75),
            "SSIM_IQR": np.nanpercentile(ssim_vals, 75) - np.nanpercentile(ssim_vals, 25),
        })

    df = pd.DataFrame(rows)

    outpath = os.path.join(folder, "plots", "metric_summary.csv")
    df.to_csv(outpath, index=False)

    print(f"✅ Metric summary saved to {outpath}")


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
            output = output.loc[ground_truth.index]
            metrics[parent_folder] = get_metrics(output, fpath)
        else:
            print(f"Shape mismatch for {parent_folder}: got {output.shape}, whereas ground truth has {ground_truth.shape}")

    plot_metrics(metrics)
    print("✅ Plots saved in 'plots/' folder.")


if __name__ == "__main__":
    main()
