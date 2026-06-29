import argparse
import os
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


METHOD_ORDER = [
    "SNMF",
    "NMF",
    "SMART",
    "CARD",
    "BayesTME",
    "SpiceMix",
    "STdeconvolve",
    "starfysh",
    "RETROFIT",
]

SPOTWISE_METRICS = [
    {
        "name": "RMSE",
        "mean": "RMSE_mean",
        "q1": "RMSE_Q1",
        "q3": "RMSE_Q3",
        "ylabel": "RMSE",
        "filename": "ngenes_rmse_comparison",
    },
    {
        "name": "JSD",
        "mean": "JSD_mean",
        "q1": "JSD_Q1",
        "q3": "JSD_Q3",
        "ylabel": "Jensen-Shannon divergence",
        "filename": "ngenes_jsd_comparison",
    },
    {
        "name": "PCC",
        "mean": "PCC_mean",
        "q1": "PCC_Q1",
        "q3": "PCC_Q3",
        "ylabel": "PCC",
        "filename": "ngenes_pcc_comparison",
    },
    {
        "name": "SSIM",
        "mean": "SSIM_mean",
        "q1": "SSIM_Q1",
        "q3": "SSIM_Q3",
        "ylabel": "SSIM",
        "filename": "ngenes_ssim_comparison",
    },
]

RUN_METRICS = [
    {
        "name": "Time",
        "value": "Time",
        "ylabel": "Time (s)",
        "filename": "ngenes_time_comparison",
        "log_scale": True,
    },
    {
        "name": "Memory",
        "value": "Memory",
        "ylabel": "Memory (MB)",
        "filename": "ngenes_memory_comparison",
        "log_scale": False,
    },
]


sns.set_theme(style="whitegrid")
sns.set_context("talk", font_scale=1.1, rc={
    "axes.titlesize": 18,
    "axes.labelsize": 15,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "axes.linewidth": 1.2,
    "grid.linewidth": 0.6,
})


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate completed ngenes ablation metric summaries and plot "
            "metrics against the number of highly variable genes."
        )
    )
    parser.add_argument(
        "results_path",
        help="Path containing v<number> gene-count result folders.",
    )
    return parser.parse_args()


def gene_count_from_dir(path):
    match = re.fullmatch(r"v(\d+)", path.name)
    if match is None:
        return None
    return int(match.group(1))


def discover_summary_files(results_path):
    summary_files = []
    skipped = []

    for child in sorted(results_path.iterdir(), key=lambda p: p.name):
        if not child.is_dir():
            continue

        ngenes = gene_count_from_dir(child)
        if ngenes is None:
            continue

        summary_path = child / "plots" / "metric_summary.csv"
        if summary_path.exists():
            summary_files.append((ngenes, summary_path))
        else:
            skipped.append((ngenes, summary_path))

    return sorted(summary_files), sorted(skipped)


def load_metric_summaries(summary_files):
    frames = []

    for ngenes, summary_path in summary_files:
        df = pd.read_csv(summary_path)
        df.insert(0, "ngenes", ngenes)
        df.insert(1, "source_summary", str(summary_path))
        frames.append(df)

    if not frames:
        return pd.DataFrame()

    return pd.concat(frames, ignore_index=True)


def ordered_methods(df):
    present = set(df["Method"].dropna().unique())
    ordered = [method for method in METHOD_ORDER if method in present]
    ordered.extend(sorted(present.difference(METHOD_ORDER)))
    return ordered


def method_color_map(methods):
    colors = sns.color_palette("Spectral", n_colors=len(methods))
    return dict(zip(methods, colors))


def validate_columns(df, required_columns):
    missing = [col for col in required_columns if col not in df.columns]
    if missing:
        raise ValueError(
            "Missing expected columns in metric summaries: "
            + ", ".join(sorted(missing))
        )


def prepare_axis(ax, ngenes_values):
    ax.set_xscale("log", base=2)
    ax.set_xticks(ngenes_values)
    ax.set_xticklabels([str(v) for v in ngenes_values])
    ax.set_xlabel("Number of highly variable genes")
    ax.grid(True, axis="y", linewidth=0.6)
    ax.grid(True, axis="x", linewidth=0.3, alpha=0.4)


def draw_legend(ax, methods, colors):
    handles = [
        plt.Line2D(
            [0],
            [0],
            color=colors[method],
            marker="o",
            linewidth=2,
            markersize=6,
            label=method,
        )
        for method in methods
    ]

    ax.legend(
        handles=handles,
        title="Method",
        title_fontsize=13,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=True,
    )


def save_current_figure(outdir, filename):
    plt.tight_layout()
    plt.savefig(outdir / f"{filename}.png", dpi=300)
    plt.savefig(outdir / f"{filename}.pdf", dpi=300)
    plt.close()


def plot_spotwise_metric(df, metric, methods, colors, ngenes_values, outdir):
    fig, ax = plt.subplots(figsize=(11, 6))

    for method in methods:
        method_df = (
            df[df["Method"] == method]
            .sort_values("ngenes")
            .dropna(subset=[metric["mean"]])
        )
        if method_df.empty:
            continue

        x = method_df["ngenes"].to_numpy(dtype=float)
        y = method_df[metric["mean"]].to_numpy(dtype=float)
        q1 = method_df[metric["q1"]].to_numpy(dtype=float)
        q3 = method_df[metric["q3"]].to_numpy(dtype=float)
        lower = np.clip(y - q1, a_min=0, a_max=None)
        upper = np.clip(q3 - y, a_min=0, a_max=None)

        ax.errorbar(
            x,
            y,
            yerr=np.vstack([lower, upper]),
            fmt="o-",
            color=colors[method],
            ecolor=colors[method],
            elinewidth=1.3,
            capsize=4,
            markersize=5.5,
            linewidth=2,
            label=method,
        )

    prepare_axis(ax, ngenes_values)
    ax.set_ylabel(metric["ylabel"])
    ax.set_title(f"{metric['name']} by Number of Highly Variable Genes", pad=15, weight="bold")
    draw_legend(ax, methods, colors)
    sns.despine()
    save_current_figure(outdir, metric["filename"])


def plot_run_metric(df, metric, methods, colors, ngenes_values, outdir):
    fig, ax = plt.subplots(figsize=(11, 6))

    for method in methods:
        method_df = (
            df[df["Method"] == method]
            .sort_values("ngenes")
            .dropna(subset=[metric["value"]])
        )
        if method_df.empty:
            continue

        ax.plot(
            method_df["ngenes"],
            method_df[metric["value"]],
            "o-",
            color=colors[method],
            markersize=5.5,
            linewidth=2,
            label=method,
        )

    prepare_axis(ax, ngenes_values)
    if metric["log_scale"]:
        ax.set_yscale("log")
    ax.set_ylabel(metric["ylabel"])
    ax.set_title(f"{metric['name']} by Number of Highly Variable Genes", pad=15, weight="bold")
    draw_legend(ax, methods, colors)
    sns.despine()
    save_current_figure(outdir, metric["filename"])


def plot_all_metrics(df, outdir):
    required_columns = ["ngenes", "Method"]
    for metric in SPOTWISE_METRICS:
        required_columns.extend([metric["mean"], metric["q1"], metric["q3"]])
    for metric in RUN_METRICS:
        required_columns.append(metric["value"])
    validate_columns(df, required_columns)

    methods = ordered_methods(df)
    colors = method_color_map(methods)
    ngenes_values = sorted(df["ngenes"].dropna().unique())

    for metric in SPOTWISE_METRICS:
        plot_spotwise_metric(df, metric, methods, colors, ngenes_values, outdir)

    for metric in RUN_METRICS:
        plot_run_metric(df, metric, methods, colors, ngenes_values, outdir)


def main():
    args = parse_args()
    results_path = Path(args.results_path)

    if not results_path.exists():
        print(f"Results path does not exist: {results_path}", file=sys.stderr)
        sys.exit(1)

    if not results_path.is_dir():
        print(f"Results path is not a directory: {results_path}", file=sys.stderr)
        sys.exit(1)

    summary_files, skipped = discover_summary_files(results_path)
    for ngenes, summary_path in skipped:
        print(f"Warning: skipping v{ngenes}; missing {summary_path}")

    if not summary_files:
        print(
            f"No completed ngenes summaries found under {results_path}",
            file=sys.stderr,
        )
        sys.exit(1)

    combined = load_metric_summaries(summary_files)
    methods = ordered_methods(combined)
    combined["_method_order"] = pd.Categorical(
        combined["Method"],
        categories=methods,
        ordered=True,
    )
    combined = (
        combined.sort_values(["ngenes", "_method_order", "Method"])
        .drop(columns="_method_order")
    )

    outdir = results_path / "plots"
    os.makedirs(outdir, exist_ok=True)

    summary_out = outdir / "ngenes_metric_summary.csv"
    combined.to_csv(summary_out, index=False)

    plot_all_metrics(combined, outdir)

    print(f"Combined metric summary saved to {summary_out}")
    print(f"Plots saved to {outdir}")


if __name__ == "__main__":
    main()
