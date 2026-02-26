
import os
import sys

import pandas as pd
import scanpy as sc

import matplotlib.pyplot as plt
import seaborn as sns

# === Configuration ===
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
    "legend.fontsize": 8,
    "axes.linewidth": 1.2,
    "grid.linewidth": 0.6
})

ground_truth = None
folder = None
hungarian = None


# === Helper functions ===
def get_csv_files():
    global folder, hungarian
    csv_files = dict()
    for root, _, files in os.walk(folder):
        for f in files:
            if f.lower().endswith(".csv") and "gene_sig" not in f and (('hungarian' in f) == hungarian):
                csv_files[os.path.basename(root)] = os.path.join(root, f)
    return csv_files

def process_proportions(proportions):
    proportions.obs[["array_row", "array_col"]] = proportions.obs_names.str.split("x").tolist()
    proportions.obs[["array_row", "array_col"]] = proportions.obs[["array_row", "array_col"]].astype(float)
    proportions.obsm['spatial'] = proportions.obs[["array_row", "array_col"]].values
    return proportions

def get_vlimits(csv_files):
    global ground_truth
    all_df = pd.DataFrame(ground_truth.X, columns=ground_truth.var_names.tolist())
    for file in csv_files.values():
        all_df = pd.concat([all_df, pd.read_csv(file, index_col=0)], axis=0)
    return {ct: (all_df[ct].min(), all_df[ct].max()) for ct in all_df.columns}

def plot_proportions(csv_files, normalize=False):
    """Plot ground truth + method proportions with shared scale, ordered rows, and labeled columns."""
    global folder, ground_truth

    os.makedirs(os.path.join(folder, "plots"), exist_ok=True)

    # Filter and order available methods by METHOD_ORDER
    available_methods = [m for m in METHOD_ORDER if m in csv_files]
    n_methods = len(available_methods)
    n_colors = len(ground_truth.var_names)
    total_panels = n_methods + 1  # +1 for ground truth

    # Layout: one row per method (plus ground truth), one column per cell type
    nrow = total_panels
    ncol = n_colors

    _, axes = plt.subplots(
        nrow,
        ncol,
        figsize=(3.0 * ncol if normalize else 4 * ncol, 3.0 * nrow),
        sharex=True,
        sharey=True,
        gridspec_kw=dict(wspace=0.05, hspace=0.05)
    )

    # Normalize axes shape
    if nrow == 1:
        axes = axes[None, :]
    if ncol == 1:
        axes = axes[:, None]

    # Shared spatial extent from ground truth
    ground_truth = process_proportions(ground_truth)
    spatial_coords = ground_truth.obsm["spatial"]
    xlim = (spatial_coords[:, 0].min(), spatial_coords[:, 0].max())
    ylim = (spatial_coords[:, 1].min(), spatial_coords[:, 1].max())

    limits = get_vlimits(csv_files)

    # === Top column titles ===
    for j, color in enumerate(ground_truth.var_names):
        axes[0, j].set_title(color, fontsize=15, pad=8)

    # === Row 0: Ground Truth ===
    for j, color in enumerate(ground_truth.var_names):
        ax = axes[0, j]
        sc.pl.spatial(ground_truth, color=color, spot_size=1.0, ax=ax, show=False,
                          vmin=limits[color][0] if normalize else None,
                          vmax=limits[color][1] if normalize else None, 
                          colorbar_loc=None if normalize else 'right')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticks([])
        ax.set_yticks([])
        if j == 0:
            ax.set_ylabel("Ground Truth", fontsize=13, rotation=0, labelpad=40, va="center")

        if not normalize:
            cbar_ax = ax.figure.axes[-1]
            cbar_ax.tick_params(
                labelsize=8,     
                width=0.5,       # thinner tick lines (default ~1–1.5)
                length=3         # shorter ticks (optional)
            )

    # === Remaining rows: methods ===
    for i, method in enumerate(available_methods, start=1):
        proportions = process_proportions(sc.read_csv(csv_files[method]))
        for j, color in enumerate(proportions.var_names):
            ax = axes[i, j]
            sc.pl.spatial(proportions, color=color, spot_size=1.0, ax=ax, show=False,
                          vmin=limits[color][0] if normalize else None,
                          vmax=limits[color][1] if normalize else None, 
                          colorbar_loc=None if normalize else 'right')
            ax.set_title('')
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.set_xticks([])
            ax.set_yticks([])
            if j == 0:
                ax.set_ylabel(method, fontsize=13, rotation=0, labelpad=40, va="center")

            if not normalize:
                cbar_ax = ax.figure.axes[-1]
                cbar_ax.tick_params(
                    labelsize=8,     
                    width=0.5,       # thinner tick lines (default ~1–1.5)
                    length=3         # shorter ticks (optional)
                )

    # plt.tight_layout(pad=0.8, w_pad=0.2, h_pad=0.2)
    plt.savefig(os.path.join(folder, "plots", f"proportions{'_norm' if normalize else ''}.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(folder, "plots", f"proportions{'_norm' if normalize else ''}.pdf"), dpi=300, bbox_inches="tight")
    plt.close()
    
def main():
    global ground_truth, folder, hungarian

    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} <folder_path> <ground_truth_path> <hungarian:true/false>")
        sys.exit(1)

    folder = sys.argv[1]
    ground_truth = process_proportions(sc.read_csv(sys.argv[2]))
    hungarian = sys.argv[3].lower() == "true"

    csv_files = get_csv_files()
    if not csv_files:
        print("No CSV files found.")
        sys.exit(0)

    plot_proportions(csv_files)
    # plot_proportions(csv_files, normalize=True)
    print("✅ Plots saved in 'plots/' folder.")


if __name__ == "__main__":
    main()
