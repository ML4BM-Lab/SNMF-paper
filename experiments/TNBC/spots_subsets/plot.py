import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# === Configuration ===
MEM_CONVERSION = {'T': 2, 'G': 1, 'M': 0, 'K': -1}

METHOD_ORDER = [
    "SMART", "STdeconvolve", "starfysh", "CARD", "RETROFIT",
    "SNMF", "SNMF_v2"
]

# Global seaborn style and consistent color palette
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
METHOD_COLOR_MAP["CARD"] = "#FFC84D"

# === Helper function ===
def parse_sacct_output(sacct_file):
    mem, time = None, None
    if os.path.exists(sacct_file):
        with open(sacct_file, 'r') as f:
            for line in f.readlines():
                if "batch" in line and "COMPLETED" in line:
                    mem, time = line.split()[2:4]
        if mem is not None:
            mem = int(mem[:-1]) * 1024 ** MEM_CONVERSION[mem[-1]]  # MB
        if time is not None:
            h, m, s = map(int, time.split(":"))
            time = h * 3600 + m * 60 + s
    return mem, time


def plot_scaling(x_values, data_dict, ylabel, title, filename, results_path, log_scale=True):
    plt.figure(figsize=(9.5, 5.5))

    methods_present = [m for m in METHOD_ORDER if any(m in d for d in data_dict.values())]
    for method in methods_present:
        y = [data_dict[s].get(method, np.nan) for s in x_values]
        plt.plot(
            x_values, y,
            marker="o", markersize=6, linewidth=2,
            label=method, color=METHOD_COLOR_MAP[method]
        )

    plt.xlabel("Number of spots", fontsize=15)
    plt.ylabel(ylabel, fontsize=15)
    plt.title(title, fontsize=18, pad=15, weight='bold')
    if log_scale:
        plt.yscale("log")

    sns.despine()
    plt.grid(True, linestyle='--', alpha=0.5)

    # Move legend outside to the right
    legend_handles = [
        plt.Line2D([0], [0], color=METHOD_COLOR_MAP[m], lw=3, marker='o', label=m)
        for m in methods_present
    ]
    plt.legend(
        handles=legend_handles,
        title="Method", title_fontsize=13,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=True
    )

    plt.tight_layout(rect=[0, 0, 0.83, 1])
    plt.savefig(os.path.join(results_path, "plots", filename), dpi=300)
    plt.close()


# === Main ===
if __name__ == "__main__":
    results_path = sys.argv[1]

    times = dict()
    mems = dict()

    # --- Collect data ---
    for root, _, files in os.walk(results_path):
        for f in files:
            if f == 'sacct.log':
                subset, method = root.split("/")[-2:]
                subset = int(subset.split("_")[0][:-5])  # adjust if needed for folder naming
                if subset not in times:
                    times[subset] = dict()
                    mems[subset] = dict()
                mems[subset][method], times[subset][method] = parse_sacct_output(os.path.join(root, f))

    subsets_sorted = sorted(times.keys())
    os.makedirs(os.path.join(results_path, "plots"), exist_ok=True)

    # --- Plot scaling metrics ---
    plot_scaling(
        subsets_sorted, times,
        ylabel="Time (s)",
        title="Runtime by Subset",
        filename="time_scaling.png",
        results_path=results_path,
        log_scale=True
    )

    plot_scaling(
        subsets_sorted, mems,
        ylabel="Memory (MB)",
        title="Memory Usage by Subset",
        filename="mem_scaling.png",
        results_path=results_path,
        log_scale=True
    )

    print("✅ Scaling plots saved: time_scaling.png, mem_scaling.png")
