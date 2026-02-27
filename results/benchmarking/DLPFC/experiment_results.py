
import sys
import os

results_path = sys.argv[1]

ari_results = {}

for (dirpath, dirnames, filenames) in os.walk(results_path):
    for file in filenames:
        if file == "ari.txt":
            sample = dirpath.split("/")[-1]
            ari_results[sample] = {}
            with open(os.path.join(dirpath, file), "r") as f:
                for i,line in enumerate(f.readlines()):
                    if i == 0:
                        continue
                    method, ari = line.split(":")
                    ari_results[sample][method] = float(ari[1:])


# Plot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
# -------------------------
# Global publication style
# -------------------------
mpl.rcParams.update({
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "font.size": 10,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "pdf.fonttype": 42,   # editable text in Illustrator
    "ps.fonttype": 42,
})

# Sort for consistent ordering
samples = sorted(ari_results.keys())
methods = [
    "SNMF",
    "STdeconvolve",
    "SpiceMix",
    "CARD",
    "SMART",
    "RETROFIT",
    "starfysh",
    "BayesTME"
]

n_samples = len(samples)
n_methods = len(methods)

# Build data matrix (rows = methods, cols = samples)
data = np.zeros((n_methods, n_samples))

for i, method in enumerate(methods):
    for j, sample in enumerate(samples):
        data[i, j] = ari_results[sample].get(method, np.nan)

# X positions
x = np.arange(n_samples)

# Better spacing
total_width = 0.7
bar_width = total_width / n_methods

# -------------------------
# Color palette
# -------------------------
colors = sns.color_palette("Spectral", n_colors=n_methods)

# -------------------------
# Create figure
# -------------------------
fig, ax = plt.subplots(figsize=(12, 5))

for i, method in enumerate(methods):
    ax.bar(
        x - total_width/2 + i * bar_width + bar_width/2,
        data[i],
        width=bar_width,
        label=method,
        color=colors[i],
        edgecolor="black",
        linewidth=0.4
    )

# Axis labels
ax.set_ylabel("Adjusted Rand Index (ARI)", fontsize=11)
# ax.set_xlabel("Sample", fontsize=11)

# X ticks centered
ax.set_xticks(x)
ax.set_xticklabels(samples)

# Subtle horizontal grid
ax.yaxis.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)
ax.set_axisbelow(True)

# Legend outside
ax.legend(
    frameon=False,
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    title="Method"
)

plt.tight_layout()

# -------------------------
# Save outputs
# -------------------------
plot_dir = os.path.join(results_path, "plots")
os.makedirs(plot_dir, exist_ok=True)

fig.savefig(os.path.join(plot_dir, "ari.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(plot_dir, "ari.svg"), bbox_inches="tight")
fig.savefig(os.path.join(plot_dir, "ari.png"), bbox_inches="tight")

plt.close(fig)

# -------------------------
# Create LaTeX table (bold best per sample)
# -------------------------

latex_path = os.path.join(results_path, "ari_table.tex")

with open(latex_path, "w") as f:
    f.write("\\begin{table*}[!t]\n")
    f.write("\\centering\n")
    f.write("\\caption{Adjusted Rand Index (ARI) for all methods across samples.}\n")
    f.write("\\label{tab:ari_results}\n")

    col_format = "l|" + "c" * n_methods
    f.write(f"\\begin{{tabular}}{{{col_format}}}\n")
    f.write("\\hline\n")

    # Header
    header = "\\textbf{Sample} & " + " & ".join([f"\\textbf{{{m}}}" for m in methods]) + " \\\\\n"
    f.write(header)
    f.write("\\hline\n")

    # Rows
    for j, sample in enumerate(samples):
        row_values = []
        
        # Extract row (all methods for this sample)
        row_data = data[:, j]
        
        # Compute best ARI ignoring NaNs
        if np.all(np.isnan(row_data)):
            best_value = None
        else:
            best_value = np.nanmax(row_data)

        for i in range(n_methods):
            value = row_data[i]

            if np.isnan(value):
                row_values.append("--")
            else:
                formatted = f"{value:.3f}"
                
                # Bold if best (handle floating point safely)
                if best_value is not None and np.isclose(value, best_value):
                    formatted = f"\\textbf{{{formatted}}}"
                
                row_values.append(formatted)

        row = f"\\textbf{{{sample}}}" + " & " + " & ".join(row_values) + " \\\\\n"
        f.write(row)

    f.write("\\hline\n")
    f.write("\\end{tabular}\n")
    f.write("\\end{table*}\n")