#!/usr/bin/env python3
import argparse
import math
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import friedmanchisquare, rankdata, studentized_range


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

BENCHMARK_DATASETS = ["TNBC", "PDAC", "Xenium"]

METRICS = {
    "RMSE": {"column": "RMSE_mean", "direction": "min", "ylabel": "RMSE"},
    "JSD": {"column": "JSD_mean", "direction": "min", "ylabel": "JSD"},
    "PCC": {"column": "PCC_mean", "direction": "max", "ylabel": "PCC"},
    "SSIM": {"column": "SSIM_mean", "direction": "max", "ylabel": "SSIM"},
}


sns.set_theme(style="whitegrid")
sns.set_context(
    "talk",
    font_scale=1.0,
    rc={
        "axes.titlesize": 17,
        "axes.labelsize": 14,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 11,
        "axes.linewidth": 1.1,
        "grid.linewidth": 0.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    },
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run Friedman-Nemenyi benchmark analysis over the existing "
            "metric_summary.csv files and summarize DLPFC ARI mean ranks."
        )
    )
    parser.add_argument(
        "benchmarking_path",
        type=Path,
        help="Path to the benchmarking folder, e.g. results_NB/benchmarking.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Family-wise alpha used for critical-difference plots.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Optional output directory. Defaults to "
            "<benchmarking_path>/plots/friedman_nemenyi."
        ),
    )
    return parser.parse_args()


def ordered_methods(methods):
    present = set(methods)
    ordered = [method for method in METHOD_ORDER if method in present]
    ordered.extend(sorted(present.difference(METHOD_ORDER)))
    return ordered


def method_color_map(methods):
    colors = sns.color_palette("Spectral", n_colors=len(methods))
    return dict(zip(methods, colors))


def metric_slug(metric_name):
    return metric_name.lower()


def save_current_figure(outdir, stem):
    plt.tight_layout()
    plt.savefig(outdir / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.savefig(outdir / f"{stem}.pdf", dpi=300, bbox_inches="tight")
    plt.close()


def load_metric_values(benchmarking_path):
    rows = []

    for dataset in BENCHMARK_DATASETS:
        summary_path = benchmarking_path / dataset / "plots" / "metric_summary.csv"
        if not summary_path.exists():
            raise FileNotFoundError(f"Missing metric summary: {summary_path}")

        summary = pd.read_csv(summary_path)
        required = {"Method", *[spec["column"] for spec in METRICS.values()]}
        missing = sorted(required.difference(summary.columns))
        if missing:
            raise ValueError(
                f"{summary_path} is missing expected columns: {', '.join(missing)}"
            )

        for _, row in summary.iterrows():
            method = row["Method"]
            for metric_name, spec in METRICS.items():
                rows.append(
                    {
                        "Dataset": dataset,
                        "Metric": metric_name,
                        "Method": method,
                        "Value": row[spec["column"]],
                        "Direction": spec["direction"],
                        "Source": str(summary_path),
                    }
                )

    values = pd.DataFrame(rows)
    values["Value"] = pd.to_numeric(values["Value"], errors="coerce")
    return values


def complete_case_matrix(values, metric_name):
    metric_values = values[values["Metric"] == metric_name]
    matrix = metric_values.pivot_table(
        index="Dataset",
        columns="Method",
        values="Value",
        aggfunc="first",
    )
    matrix = matrix.reindex(BENCHMARK_DATASETS)

    observed_methods = ordered_methods(matrix.columns.dropna().tolist())
    included_methods = [
        method
        for method in observed_methods
        if method in matrix.columns and matrix[method].notna().all()
    ]
    dropped_methods = [method for method in observed_methods if method not in included_methods]

    if len(included_methods) < 2:
        raise ValueError(
            f"{metric_name} has fewer than two complete-case methods: "
            f"{included_methods}"
        )

    matrix = matrix[included_methods]
    if matrix.isna().any().any():
        raise ValueError(f"{metric_name} complete-case matrix still contains NaN values")

    return matrix, included_methods, dropped_methods


def rank_metric_matrix(matrix, metric_name):
    direction = METRICS[metric_name]["direction"]
    rank_rows = []
    rank_matrix = pd.DataFrame(index=matrix.index, columns=matrix.columns, dtype=float)

    for dataset, row in matrix.iterrows():
        values = row.to_numpy(dtype=float)
        oriented = values if direction == "min" else -values
        ranks = rankdata(oriented, method="average")

        for method, value, rank in zip(matrix.columns, values, ranks):
            rank_matrix.loc[dataset, method] = rank
            rank_rows.append(
                {
                    "Dataset": dataset,
                    "Metric": metric_name,
                    "Method": method,
                    "Value": value,
                    "Rank": rank,
                    "Direction": direction,
                }
            )

    return rank_matrix, pd.DataFrame(rank_rows)


def nemenyi_pvalues(average_ranks, n_blocks):
    methods = list(average_ranks.index)
    n_methods = len(methods)
    standard_error = math.sqrt(n_methods * (n_methods + 1) / (6.0 * n_blocks))
    pvalues = pd.DataFrame(1.0, index=methods, columns=methods, dtype=float)

    for i, method_i in enumerate(methods):
        for j, method_j in enumerate(methods):
            if i == j:
                continue
            q_stat = abs(average_ranks[method_i] - average_ranks[method_j]) / standard_error
            pvalues.loc[method_i, method_j] = studentized_range.sf(
                q_stat * math.sqrt(2.0),
                n_methods,
                np.inf,
            )

    return pvalues.clip(lower=0.0, upper=1.0)


def critical_difference(n_methods, n_blocks, alpha):
    q_alpha = studentized_range.ppf(1.0 - alpha, n_methods, np.inf) / math.sqrt(2.0)
    return q_alpha * math.sqrt(n_methods * (n_methods + 1) / (6.0 * n_blocks))


def friedman_nemenyi(values, alpha):
    rank_frames = []
    average_rank_frames = []
    friedman_rows = []
    nemenyi_results = {}

    for metric_name in METRICS:
        matrix, included_methods, dropped_methods = complete_case_matrix(values, metric_name)
        rank_matrix, rank_long = rank_metric_matrix(matrix, metric_name)
        rank_frames.append(rank_long)

        statistic, p_value = friedmanchisquare(
            *[matrix[method].to_numpy(dtype=float) for method in included_methods]
        )

        average_ranks = rank_matrix.mean(axis=0).sort_values()
        average_values = matrix.mean(axis=0).reindex(average_ranks.index)
        rank_sd = rank_matrix.std(axis=0, ddof=0).reindex(average_ranks.index)
        n_blocks = matrix.shape[0]
        n_methods = matrix.shape[1]
        cd = critical_difference(n_methods, n_blocks, alpha)

        for order, method in enumerate(average_ranks.index, start=1):
            average_rank_frames.append(
                {
                    "Metric": metric_name,
                    "Method": method,
                    "Average_rank": average_ranks[method],
                    "Rank_sd": rank_sd[method],
                    "Mean_value": average_values[method],
                    "Order": order,
                    "N_datasets": n_blocks,
                    "Direction": METRICS[metric_name]["direction"],
                }
            )

        friedman_rows.append(
            {
                "Metric": metric_name,
                "Statistic": statistic,
                "P_value": p_value,
                "N_datasets": n_blocks,
                "N_methods": n_methods,
                "Alpha": alpha,
                "Critical_difference": cd,
                "Included_methods": ";".join(included_methods),
                "Dropped_methods": ";".join(dropped_methods),
                "Datasets": ";".join(matrix.index.tolist()),
            }
        )

        nemenyi_results[metric_name] = nemenyi_pvalues(average_ranks, n_blocks)

    ranks = pd.concat(rank_frames, ignore_index=True)
    average_ranks = pd.DataFrame(average_rank_frames)
    friedman = pd.DataFrame(friedman_rows)
    return ranks, average_ranks, friedman, nemenyi_results


def build_collapsed_block_values(metric_values, dlpfc_ranks):
    rows = []

    for _, row in metric_values.iterrows():
        rows.append(
            {
                "Block": f"{row['Dataset']}_{row['Metric']}",
                "Dataset": row["Dataset"],
                "Metric": row["Metric"],
                "Method": row["Method"],
                "Value": row["Value"],
                "Direction": row["Direction"],
                "Source": row["Source"],
            }
        )

    for _, row in dlpfc_ranks.iterrows():
        rows.append(
            {
                "Block": "DLPFC_ARI",
                "Dataset": "DLPFC",
                "Metric": "ARI",
                "Method": row["Method"],
                "Value": row["Value"],
                "Direction": "max",
                "Source": "DLPFC/results/ari_table.tex:Mean",
            }
        )

    block_values = pd.DataFrame(rows)
    block_values["Value"] = pd.to_numeric(block_values["Value"], errors="coerce")
    return block_values


def collapsed_complete_case_matrix(block_values):
    block_order = [
        f"{dataset}_{metric}"
        for dataset in BENCHMARK_DATASETS
        for metric in METRICS
    ]
    block_order.append("DLPFC_ARI")

    block_info = (
        block_values[["Block", "Dataset", "Metric", "Direction"]]
        .drop_duplicates()
        .set_index("Block")
        .reindex(block_order)
    )

    if block_info.isna().any().any():
        missing = block_info[block_info.isna().any(axis=1)].index.tolist()
        raise ValueError(f"Missing collapsed-test blocks: {', '.join(missing)}")

    matrix = block_values.pivot_table(
        index="Block",
        columns="Method",
        values="Value",
        aggfunc="first",
    )
    matrix = matrix.reindex(block_order)

    observed_methods = ordered_methods(matrix.columns.dropna().tolist())
    included_methods = [
        method
        for method in observed_methods
        if method in matrix.columns and matrix[method].notna().all()
    ]
    dropped_methods = [method for method in observed_methods if method not in included_methods]

    if len(included_methods) < 2:
        raise ValueError(
            "Collapsed accuracy test has fewer than two complete-case methods: "
            f"{included_methods}"
        )

    matrix = matrix[included_methods]
    if matrix.isna().any().any():
        raise ValueError("Collapsed complete-case matrix still contains NaN values")

    return matrix, block_info, included_methods, dropped_methods


def rank_collapsed_matrix(matrix, block_info):
    rank_rows = []
    rank_matrix = pd.DataFrame(index=matrix.index, columns=matrix.columns, dtype=float)

    for block, row in matrix.iterrows():
        direction = block_info.loc[block, "Direction"]
        values = row.to_numpy(dtype=float)
        oriented = values if direction == "min" else -values
        ranks = rankdata(oriented, method="average")

        for method, value, rank in zip(matrix.columns, values, ranks):
            rank_matrix.loc[block, method] = rank
            rank_rows.append(
                {
                    "Block": block,
                    "Dataset": block_info.loc[block, "Dataset"],
                    "Metric": block_info.loc[block, "Metric"],
                    "Method": method,
                    "Value": value,
                    "Rank": rank,
                    "Direction": direction,
                }
            )

    return rank_matrix, pd.DataFrame(rank_rows)


def collapsed_friedman_nemenyi(block_values, alpha):
    matrix, block_info, included_methods, dropped_methods = collapsed_complete_case_matrix(
        block_values
    )
    rank_matrix, rank_long = rank_collapsed_matrix(matrix, block_info)

    statistic, p_value = friedmanchisquare(
        *[rank_matrix[method].to_numpy(dtype=float) for method in included_methods]
    )

    average_ranks = rank_matrix.mean(axis=0).sort_values()
    rank_sd = rank_matrix.std(axis=0, ddof=0).reindex(average_ranks.index)
    n_blocks = rank_matrix.shape[0]
    n_methods = rank_matrix.shape[1]
    cd = critical_difference(n_methods, n_blocks, alpha)
    pvalues = nemenyi_pvalues(average_ranks, n_blocks)

    average_rows = []
    for order, method in enumerate(average_ranks.index, start=1):
        average_rows.append(
            {
                "Test": "collapsed_accuracy_rank",
                "Method": method,
                "Average_rank": average_ranks[method],
                "Rank_sd": rank_sd[method],
                "Order": order,
                "N_blocks": n_blocks,
            }
        )

    friedman = pd.DataFrame(
        [
            {
                "Test": "collapsed_accuracy_rank",
                "Statistic": statistic,
                "P_value": p_value,
                "N_blocks": n_blocks,
                "N_methods": n_methods,
                "Alpha": alpha,
                "Critical_difference": cd,
                "Included_methods": ";".join(included_methods),
                "Dropped_methods": ";".join(dropped_methods),
                "Blocks": ";".join(rank_matrix.index.tolist()),
            }
        ]
    )

    return rank_long, pd.DataFrame(average_rows), friedman, pvalues


def clean_latex_cell(cell):
    cell = cell.strip()
    cell = re.sub(r"\\textbf\{([^{}]+)\}", r"\1", cell)
    cell = cell.replace("\\\\", "")
    cell = cell.replace("$", "")
    return cell.strip()


def parse_dlpfc_ari_mean(benchmarking_path):
    ari_table = benchmarking_path / "DLPFC" / "results" / "ari_table.tex"
    if not ari_table.exists():
        raise FileNotFoundError(f"Missing DLPFC ARI table: {ari_table}")

    lines = ari_table.read_text().splitlines()
    header_line = next(
        (
            line
            for line in lines
            if "\\textbf{Sample}" in line and "&" in line and "\\textbf{Mean}" not in line
        ),
        None,
    )
    mean_line = next((line for line in lines if line.strip().startswith("\\textbf{Mean}")), None)

    if header_line is None or mean_line is None:
        raise ValueError(f"Could not parse header and Mean rows from {ari_table}")

    methods = [clean_latex_cell(cell) for cell in header_line.split("&")[1:]]
    values = []
    for cell in mean_line.split("&")[1:]:
        clean = clean_latex_cell(cell)
        values.append(np.nan if clean == "--" else float(clean))

    if len(methods) != len(values):
        raise ValueError(
            f"Parsed {len(methods)} DLPFC methods but {len(values)} mean values"
        )

    df = pd.DataFrame(
        {
            "Dataset": "DLPFC",
            "Metric": "ARI",
            "Method": methods,
            "Value": values,
            "Direction": "max",
        }
    ).dropna(subset=["Value"])

    df["Rank"] = rankdata(-df["Value"].to_numpy(dtype=float), method="average")
    df = df.sort_values(["Rank", "Method"]).reset_index(drop=True)
    return df


def plot_average_ranks(average_ranks, outdir, methods):
    colors = method_color_map(methods)

    for metric_name in METRICS:
        metric_df = average_ranks[average_ranks["Metric"] == metric_name].sort_values(
            "Average_rank"
        )
        if metric_df.empty:
            continue

        fig, ax = plt.subplots(figsize=(9.5, 5.5))
        bar_colors = [colors[method] for method in metric_df["Method"]]
        ax.barh(metric_df["Method"], metric_df["Average_rank"], color=bar_colors, edgecolor="black")
        ax.invert_yaxis()
        ax.set_xlabel("Average rank (lower is better)")
        ax.set_ylabel("")
        ax.set_title(f"{metric_name}: average Friedman ranks")
        ax.set_xlim(0, max(metric_df["Average_rank"].max() + 0.75, 2.0))
        ax.grid(True, axis="x", linewidth=0.6)
        ax.grid(False, axis="y")

        for y_pos, (_, row) in enumerate(metric_df.iterrows()):
            ax.text(
                row["Average_rank"] + 0.08,
                y_pos,
                f"{row['Average_rank']:.2f}",
                va="center",
                ha="left",
                fontsize=10,
            )

        sns.despine()
        save_current_figure(outdir, f"average_ranks_{metric_name.lower()}")


def plot_nemenyi_heatmaps(nemenyi_results, outdir):
    for metric_name, pvalues in nemenyi_results.items():
        fig, ax = plt.subplots(figsize=(8.2, 7.0))
        sns.heatmap(
            pvalues,
            ax=ax,
            vmin=0,
            vmax=1,
            cmap="rocket_r",
            square=True,
            annot=True,
            fmt=".2g",
            linewidths=0.4,
            linecolor="white",
            annot_kws={"fontsize": 8},
            cbar_kws={"label": "Nemenyi p-value"},
        )
        ax.set_title(f"{metric_name}: pairwise Nemenyi p-values")
        ax.set_xlabel("")
        ax.set_ylabel("")
        plt.xticks(rotation=45, ha="right")
        plt.yticks(rotation=0)
        ax.tick_params(axis="both", labelsize=10)
        ax.figure.axes[-1].tick_params(labelsize=10)
        ax.figure.axes[-1].yaxis.label.set_size(11)
        save_current_figure(outdir, f"nemenyi_pvalues_{metric_name.lower()}_heatmap")


def plot_critical_difference(average_ranks, friedman_results, outdir, methods):
    colors = method_color_map(methods)

    for metric_name in METRICS:
        metric_df = average_ranks[average_ranks["Metric"] == metric_name].sort_values(
            "Average_rank"
        )
        if metric_df.empty:
            continue
        friedman_row = friedman_results[friedman_results["Metric"] == metric_name].iloc[0]
        cd = friedman_row["Critical_difference"]
        n_methods = int(friedman_row["N_methods"])

        fig, ax = plt.subplots(figsize=(10, 4.8))
        y_positions = np.arange(len(metric_df))

        for y_pos, (_, row) in zip(y_positions, metric_df.iterrows()):
            method = row["Method"]
            rank = row["Average_rank"]
            ax.scatter(
                rank,
                y_pos,
                s=95,
                color=colors[method],
                edgecolor="black",
                linewidth=0.8,
                zorder=3,
            )
            ax.text(rank + 0.08, y_pos, method, va="center", ha="left", fontsize=11)

        cd_y = -0.9
        cd_start = 1.0
        cd_end = min(n_methods + 0.5, cd_start + cd)
        ax.hlines(cd_y, cd_start, cd_end, color="black", linewidth=2)
        ax.vlines([cd_start, cd_end], cd_y - 0.08, cd_y + 0.08, color="black", linewidth=2)
        ax.text(
            (cd_start + cd_end) / 2,
            cd_y,
            f"CD = {cd:.2f}",
            ha="center",
            va="center",
            fontsize=10,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5},
        )

        ax.set_xlim(0.5, n_methods + 1.0)
        ax.set_ylim(len(metric_df) - 0.4, -1.35)
        ax.set_xticks(range(1, n_methods + 1))
        ax.set_yticks([])
        ax.set_xlabel("Average rank (lower is better)")
        ax.set_title(f"{metric_name}: critical-difference rank plot")
        ax.grid(True, axis="x", linewidth=0.5, alpha=0.6)
        ax.grid(False, axis="y")
        sns.despine(left=True)
        save_current_figure(outdir, f"critical_difference_{metric_name.lower()}")


def plot_dlpfc_ari(dlpfc_ranks, outdir):
    methods = dlpfc_ranks["Method"].tolist()
    colors = method_color_map(methods)
    plot_df = dlpfc_ranks.sort_values("Rank")

    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    bar_colors = [colors[method] for method in plot_df["Method"]]
    ax.barh(plot_df["Method"], plot_df["Rank"], color=bar_colors, edgecolor="black")
    ax.invert_yaxis()
    ax.set_xlabel("ARI mean rank (lower is better)")
    ax.set_ylabel("")
    ax.set_title("DLPFC: ARI mean ranks")
    ax.set_xlim(0, max(plot_df["Rank"].max() + 0.75, 2.0))
    ax.grid(True, axis="x", linewidth=0.6)
    ax.grid(False, axis="y")

    for y_pos, (_, row) in enumerate(plot_df.iterrows()):
        ax.text(
            row["Rank"] + 0.08,
            y_pos,
            f"rank {row['Rank']:.0f}; ARI {row['Value']:.3f}",
            va="center",
            ha="left",
            fontsize=10,
        )

    sns.despine()
    save_current_figure(outdir, "dlpfc_ari_mean_ranks")


def plot_collapsed_average_ranks(average_ranks, outdir):
    plot_df = average_ranks.sort_values("Average_rank")
    methods = plot_df["Method"].tolist()
    colors = method_color_map(methods)

    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    bar_colors = [colors[method] for method in plot_df["Method"]]
    ax.barh(plot_df["Method"], plot_df["Average_rank"], color=bar_colors, edgecolor="black")
    ax.invert_yaxis()
    ax.set_xlabel("Average rank (lower is better)")
    ax.set_ylabel("")
    ax.set_title("Collapsed accuracy: average Friedman ranks")
    ax.set_xlim(0, max(plot_df["Average_rank"].max() + 0.75, 2.0))
    ax.grid(True, axis="x", linewidth=0.6)
    ax.grid(False, axis="y")

    for y_pos, (_, row) in enumerate(plot_df.iterrows()):
        ax.text(
            row["Average_rank"] + 0.08,
            y_pos,
            f"{row['Average_rank']:.2f}",
            va="center",
            ha="left",
            fontsize=10,
        )

    sns.despine()
    save_current_figure(outdir, "collapsed_average_ranks")


def plot_collapsed_heatmap(pvalues, outdir):
    fig, ax = plt.subplots(figsize=(8.2, 7.0))
    sns.heatmap(
        pvalues,
        ax=ax,
        vmin=0,
        vmax=1,
        cmap="rocket_r",
        square=True,
        annot=True,
        fmt=".2g",
        linewidths=0.4,
        linecolor="white",
        annot_kws={"fontsize": 8},
        cbar_kws={"label": "Nemenyi p-value"},
    )
    ax.set_title("Collapsed accuracy: pairwise Nemenyi p-values")
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    ax.tick_params(axis="both", labelsize=10)
    ax.figure.axes[-1].tick_params(labelsize=10)
    ax.figure.axes[-1].yaxis.label.set_size(11)
    save_current_figure(outdir, "collapsed_nemenyi_pvalues_heatmap")


def plot_collapsed_critical_difference(average_ranks, friedman_results, outdir):
    plot_df = average_ranks.sort_values("Average_rank")
    methods = plot_df["Method"].tolist()
    colors = method_color_map(methods)
    friedman_row = friedman_results.iloc[0]
    cd = friedman_row["Critical_difference"]
    n_methods = int(friedman_row["N_methods"])

    fig, ax = plt.subplots(figsize=(10, 4.8))
    y_positions = np.arange(len(plot_df))

    for y_pos, (_, row) in zip(y_positions, plot_df.iterrows()):
        method = row["Method"]
        rank = row["Average_rank"]
        ax.scatter(
            rank,
            y_pos,
            s=95,
            color=colors[method],
            edgecolor="black",
            linewidth=0.8,
            zorder=3,
        )
        ax.text(rank + 0.08, y_pos, method, va="center", ha="left", fontsize=11)

    cd_y = -0.9
    cd_start = 1.0
    cd_end = min(n_methods + 0.5, cd_start + cd)
    ax.hlines(cd_y, cd_start, cd_end, color="black", linewidth=2)
    ax.vlines([cd_start, cd_end], cd_y - 0.08, cd_y + 0.08, color="black", linewidth=2)
    ax.text(
        (cd_start + cd_end) / 2,
        cd_y,
        f"CD = {cd:.2f}",
        ha="center",
        va="center",
        fontsize=10,
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5},
    )

    ax.set_xlim(0.5, n_methods + 1.0)
    ax.set_ylim(len(plot_df) - 0.4, -1.35)
    ax.set_xticks(range(1, n_methods + 1))
    ax.set_yticks([])
    ax.set_xlabel("Average rank (lower is better)")
    ax.set_title("Collapsed accuracy: critical-difference rank plot")
    ax.grid(True, axis="x", linewidth=0.5, alpha=0.6)
    ax.grid(False, axis="y")
    sns.despine(left=True)
    save_current_figure(outdir, "collapsed_critical_difference")


def save_per_metric_outputs(
    outdir,
    metric_values,
    metric_ranks,
    average_ranks,
    friedman_results,
    nemenyi_results,
    dlpfc_ranks,
):
    outdir.mkdir(parents=True, exist_ok=True)

    metric_manifest_rows = []

    metric_values.to_csv(outdir / "metric_values_long.csv", index=False)
    metric_ranks.to_csv(outdir / "metric_ranks_long.csv", index=False)
    average_ranks.to_csv(outdir / "average_ranks.csv", index=False)
    friedman_results.to_csv(outdir / "friedman_results.csv", index=False)
    dlpfc_ranks.to_csv(outdir / "dlpfc_ari_mean_ranks.csv", index=False)

    for metric_name, pvalues in nemenyi_results.items():
        pvalues.to_csv(outdir / f"nemenyi_pvalues_{metric_name.lower()}.csv")

        metric_dir = outdir / metric_slug(metric_name)
        metric_dir.mkdir(parents=True, exist_ok=True)

        metric_values_df = metric_values[metric_values["Metric"] == metric_name]
        metric_ranks_df = metric_ranks[metric_ranks["Metric"] == metric_name]
        metric_average_df = average_ranks[average_ranks["Metric"] == metric_name]
        metric_friedman_df = friedman_results[friedman_results["Metric"] == metric_name]

        metric_values_df.to_csv(metric_dir / "metric_values.csv", index=False)
        metric_ranks_df.to_csv(metric_dir / "metric_ranks.csv", index=False)
        metric_average_df.to_csv(metric_dir / "average_ranks.csv", index=False)
        metric_friedman_df.to_csv(metric_dir / "friedman_results.csv", index=False)
        pvalues.to_csv(metric_dir / "nemenyi_pvalues.csv")

        metric_methods = ordered_methods(metric_average_df["Method"].dropna().unique().tolist())
        plot_average_ranks(metric_average_df, metric_dir, metric_methods)
        plot_nemenyi_heatmaps({metric_name: pvalues}, metric_dir)
        plot_critical_difference(
            metric_average_df,
            metric_friedman_df,
            metric_dir,
            metric_methods,
        )

        metric_manifest_rows.append(
            {
                "Metric": metric_name,
                "Description": f"Per-metric Friedman-Nemenyi test for {metric_name}.",
                "Output_dir": str(metric_dir),
            }
        )

    dlpfc_dir = outdir / "dlpfc_ari"
    dlpfc_dir.mkdir(parents=True, exist_ok=True)
    dlpfc_ranks.to_csv(dlpfc_dir / "dlpfc_ari_mean_ranks.csv", index=False)
    plot_dlpfc_ari(dlpfc_ranks, dlpfc_dir)
    metric_manifest_rows.append(
        {
            "Metric": "DLPFC_ARI",
            "Description": "DLPFC ARI descriptive mean-row ranking.",
            "Output_dir": str(dlpfc_dir),
        }
    )

    pd.DataFrame(metric_manifest_rows).to_csv(
        outdir / "per_metric_manifest.csv",
        index=False,
    )


def save_collapsed_outputs(
    outdir,
    block_values,
    rank_long,
    average_ranks,
    friedman_results,
    pvalues,
):
    outdir.mkdir(parents=True, exist_ok=True)

    block_values.to_csv(outdir / "collapsed_metric_values_long.csv", index=False)
    rank_long.to_csv(outdir / "collapsed_metric_ranks_long.csv", index=False)
    average_ranks.to_csv(outdir / "collapsed_average_ranks.csv", index=False)
    friedman_results.to_csv(outdir / "collapsed_friedman_results.csv", index=False)
    pvalues.to_csv(outdir / "collapsed_nemenyi_pvalues.csv")

    plot_collapsed_average_ranks(average_ranks, outdir)
    plot_collapsed_heatmap(pvalues, outdir)
    plot_collapsed_critical_difference(average_ranks, friedman_results, outdir)


def save_manifest(outdir, per_metric_dir, collapsed_dir):
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = pd.DataFrame(
        [
            {
                "Test": "per_metric_friedman_nemenyi",
                "Description": "Separate Friedman-Nemenyi tests for RMSE, JSD, PCC, and SSIM.",
                "Output_dir": str(per_metric_dir),
            },
            {
                "Test": "collapsed_accuracy_friedman_nemenyi",
                "Description": (
                    "Single rank-based Friedman-Nemenyi test over TNBC/PDAC/Xenium "
                    "metric blocks plus the DLPFC ARI mean block."
                ),
                "Output_dir": str(collapsed_dir),
            },
        ]
    )
    manifest.to_csv(outdir / "analysis_manifest.csv", index=False)


def main():
    args = parse_args()
    benchmarking_path = args.benchmarking_path
    outdir = args.output_dir or benchmarking_path / "plots" / "friedman_nemenyi"
    per_metric_dir = outdir / "per_metric_friedman_nemenyi"
    collapsed_dir = outdir / "collapsed_accuracy_friedman_nemenyi"

    metric_values = load_metric_values(benchmarking_path)
    metric_ranks, average_ranks, friedman_results, nemenyi_results = friedman_nemenyi(
        metric_values,
        args.alpha,
    )
    dlpfc_ranks = parse_dlpfc_ari_mean(benchmarking_path)
    collapsed_values = build_collapsed_block_values(metric_values, dlpfc_ranks)
    (
        collapsed_ranks,
        collapsed_average_ranks,
        collapsed_friedman_results,
        collapsed_nemenyi_pvalues,
    ) = collapsed_friedman_nemenyi(collapsed_values, args.alpha)

    save_per_metric_outputs(
        per_metric_dir,
        metric_values,
        metric_ranks,
        average_ranks,
        friedman_results,
        nemenyi_results,
        dlpfc_ranks,
    )
    save_collapsed_outputs(
        collapsed_dir,
        collapsed_values,
        collapsed_ranks,
        collapsed_average_ranks,
        collapsed_friedman_results,
        collapsed_nemenyi_pvalues,
    )
    save_manifest(outdir, per_metric_dir, collapsed_dir)

    print(f"Saved per-metric Friedman-Nemenyi analysis to {per_metric_dir}")
    print(f"Saved collapsed accuracy Friedman-Nemenyi analysis to {collapsed_dir}")


if __name__ == "__main__":
    main()
