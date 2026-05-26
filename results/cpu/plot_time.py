import sys
import os
import re

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
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
    return "-"


def slurm_time_to_seconds(t):
    """
    Parses Slurm-like times:
      MM:SS
      HH:MM:SS
      D-HH:MM:SS
    """
    t = t.strip()

    days = 0
    if "-" in t:
        d, t = t.split("-", 1)
        days = int(d)

    parts = [int(x) for x in t.split(":")]

    if len(parts) == 2:
        hours = 0
        minutes, seconds = parts
    elif len(parts) == 3:
        hours, minutes, seconds = parts
    else:
        raise ValueError(f"Unrecognized time format: {t}")

    return days * 86400 + hours * 3600 + minutes * 60 + seconds


def add_significance(ax, data, x_col, y_col, order, pairs, alternative="two-sided"):
    ymax = data[y_col].max()
    y_offset = ymax * 0.08
    height = ymax * 0.03
    current_y = ymax + y_offset

    for g1, g2 in pairs:
        d1 = data[data[x_col] == g1][y_col]
        d2 = data[data[x_col] == g2][y_col]

        _, p = mannwhitneyu(d1, d2, alternative=alternative)
        stars = pvalue_to_stars(p)

        x1 = order.index(g1)
        x2 = order.index(g2)

        ax.plot(
            [x1, x1, x2, x2],
            [current_y, current_y + height, current_y + height, current_y],
            lw=1.5,
            c="black",
        )

        ax.text(
            (x1 + x2) / 2,
            current_y + height,
            stars,
            ha="center",
            va="bottom",
            fontsize=14,
            fontweight="bold",
        )

        current_y += y_offset

    ax.set_ylim(top=current_y + y_offset)


sns.set_theme(style="whitegrid")
sns.set_context(
    "talk",
    font_scale=1.1,
    rc={
        "axes.titlesize": 18,
        "axes.labelsize": 15,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 12,
        "axes.linewidth": 1.2,
        "grid.linewidth": 0.6,
    },
)

results_path = sys.argv[1]
outdir = sys.argv[2] if len(sys.argv) > 2 else results_path
os.makedirs(outdir, exist_ok=True)

times = {
    "cpu": [],
    "gpu": [],
}

for dirpath, subdirs, files in os.walk(results_path):
    for filename in files:
        if filename == "sacct.log":
            pu = dirpath.split("/")[-2].lower()

            if pu not in times:
                continue

            with open(os.path.join(dirpath, filename), "r") as fh:
                for line in fh:
                    if "batch" in line and "COMPLETED" in line:
                        elapsed = line.split()[3]
                        times[pu].append(elapsed)


rows = []
for pu, vals in times.items():
    for i, t in enumerate(vals):
        seconds = slurm_time_to_seconds(t)
        rows.append(
            {
                "processing_unit": pu.upper(),
                "seed": i,
                "time_raw": t,
                "time_seconds": seconds,
                "time_minutes": seconds / 60,
                "time_hours": seconds / 3600,
            }
        )

df = pd.DataFrame(rows)

if df.empty:
    raise RuntimeError("No completed batch jobs found in sacct.log files.")

csv_path = os.path.join(outdir, "computing_times.csv")
df.to_csv(csv_path, index=False)

order = ["CPU", "GPU"]
order = [x for x in order if x in df["processing_unit"].unique()]

fig, ax = plt.subplots(figsize=(7, 6))

sns.barplot(
    data=df,
    x="processing_unit",
    y="time_minutes",
    order=order,
    errorbar="sd",
    capsize=0.15,
    width=0.6,
    edgecolor="black",
    linewidth=1.2,
    ax=ax,
)

sns.stripplot(
    data=df,
    x="processing_unit",
    y="time_minutes",
    order=order,
    color="black",
    alpha=0.7,
    size=7,
    jitter=0.12,
    ax=ax,
)

if len(order) == 2:
    add_significance(
        ax=ax,
        data=df,
        x_col="processing_unit",
        y_col="time_minutes",
        order=order,
        pairs=[("CPU", "GPU")],
        alternative="two-sided",
    )

ax.set_title("Computing time across seeds", fontweight="bold", pad=15)
ax.set_xlabel("")
ax.set_ylabel("Time / min")

sns.despine(ax=ax)
plt.tight_layout()

png_path = os.path.join(outdir, "computing_times_barplot.png")
pdf_path = os.path.join(outdir, "computing_times_barplot.pdf")

plt.savefig(png_path, dpi=300, bbox_inches="tight")
plt.savefig(pdf_path, bbox_inches="tight")
plt.close()

print(f"Saved: {csv_path}")
print(f"Saved: {png_path}")
print(f"Saved: {pdf_path}")