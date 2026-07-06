#!/usr/bin/env python3

import argparse
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
from scipy.stats import hypergeom, spearmanr


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULTS_PATH = REPO_ROOT / "results_NB" / "ablation_study" / "tau"
DEFAULT_MARKER_PATHS = {
    "TNBC": REPO_ROOT / "data" / "TNBC" / "processed" / "marker_genes.csv",
    "PDAC": REPO_ROOT
    / "data"
    / "scDesign3"
    / "PDAC"
    / "marker_genes"
    / "marker_genes.csv",
    "Xenium": REPO_ROOT
    / "data"
    / "Xenium"
    / "STHELAR"
    / "sdata_lung_s1"
    / "pseudospots"
    / "marker_genes.csv",
}
DEFAULT_FULL_MARKER_PATHS = {
    "PDAC": REPO_ROOT
    / "data"
    / "scDesign3"
    / "PDAC"
    / "marker_genes"
    / "marker_genes_full.csv",
    "Xenium": REPO_ROOT
    / "data"
    / "Xenium"
    / "STHELAR"
    / "sdata_lung_s1"
    / "pseudospots"
    / "marker_genes_full.csv",
}

FULL_MARKER_LOGFC_COLUMNS = [
    "avg_log2FC",
    "avg_logFC",
    "logfoldchange",
    "logFC",
    "log2FC",
]
FULL_MARKER_P_ADJ_COLUMNS = [
    "p_val_adj",
    "pvals_adj",
    "p_adj",
    "padj",
    "adj_p_val",
]
FULL_MARKER_P_VALUE_COLUMNS = [
    "p_val",
    "pvals",
    "p_value",
    "pvalue",
]
FULL_MARKER_SCORE_COLUMNS = ["score", "scores"]
P_ADJ_FLOOR = 1e-300


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compare tau ablation signature matrices against plain marker-gene lists."
        )
    )
    parser.add_argument(
        "--results-path",
        type=Path,
        default=DEFAULT_RESULTS_PATH,
        help=f"Tau ablation results directory. Default: {DEFAULT_RESULTS_PATH}",
    )
    parser.add_argument(
        "--dataset",
        choices=["all", *DEFAULT_MARKER_PATHS.keys()],
        default="all",
        help="Dataset to analyze. Default: all",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=50,
        help="Number of top signature genes to test per cell type. Default: 50",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="FDR threshold used to flag enriched signatures. Default: 0.05",
    )
    parser.add_argument(
        "--marker-mode",
        choices=["plain", "full", "both"],
        default="plain",
        help=(
            "Use plain marker_genes.csv, full marker_genes_full.csv DE tables, "
            "or both. Full mode is configured for PDAC and Xenium. Default: plain"
        ),
    )
    parser.add_argument(
        "--marker-p-adj",
        type=float,
        default=0.05,
        help="Adjusted p-value cutoff for --marker-mode full. Default: 0.05",
    )
    parser.add_argument(
        "--marker-logfc-min",
        type=float,
        default=0.0,
        help="Minimum logFC for --marker-mode full. Default: 0.0",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Write CSV outputs without generating plots.",
    )
    return parser.parse_args()


def normalize_label(label):
    label = str(label).strip().lower()
    label = re.sub(r"[^0-9a-z]+", "_", label)
    label = re.sub(r"_+", "_", label)
    return label.strip("_")


def read_plain_marker_genes(marker_path):
    markers = pd.read_csv(marker_path)
    missing_cols = {"gene", "cluster"} - set(markers.columns)
    if missing_cols:
        raise ValueError(
            f"{marker_path} is missing required columns: {sorted(missing_cols)}"
        )

    markers = markers[["gene", "cluster"]].dropna()
    markers["gene"] = markers["gene"].astype(str).str.strip()
    markers["cluster"] = markers["cluster"].astype(str).str.strip()
    markers = markers[(markers["gene"] != "") & (markers["cluster"] != "")]

    grouped = (
        markers.groupby("cluster")["gene"]
        .apply(lambda values: sorted(set(values)))
        .to_dict()
    )
    if not grouped:
        raise ValueError(f"No marker genes found in {marker_path}")

    normalized = normalized_cluster_map(grouped, marker_path)
    marker_details = {
        cluster: {gene: {} for gene in genes}
        for cluster, genes in grouped.items()
    }

    return grouped, normalized, marker_details


def read_full_marker_genes(marker_path, p_adj_cutoff, logfc_min):
    markers = pd.read_csv(marker_path)
    missing_cols = {"gene", "cluster"} - set(markers.columns)
    if missing_cols:
        raise ValueError(
            f"{marker_path} is missing required columns: {sorted(missing_cols)}"
        )

    logfc_col = find_column(markers, FULL_MARKER_LOGFC_COLUMNS, marker_path)
    p_adj_col = find_column(markers, FULL_MARKER_P_ADJ_COLUMNS, marker_path)
    p_value_col = find_column(
        markers,
        FULL_MARKER_P_VALUE_COLUMNS,
        marker_path,
        required=False,
    )
    score_col = find_column(
        markers,
        FULL_MARKER_SCORE_COLUMNS,
        marker_path,
        required=False,
    )

    keep_cols = ["gene", "cluster", logfc_col, p_adj_col]
    optional_cols = [col for col in [p_value_col, score_col] if col is not None]
    markers = markers[keep_cols + optional_cols].copy()
    markers = markers.rename(
        columns={
            logfc_col: "marker_logfc",
            p_adj_col: "marker_p_adj",
            p_value_col: "marker_p_value" if p_value_col else p_value_col,
            score_col: "marker_score" if score_col else score_col,
        }
    )

    for col in ["marker_p_value", "marker_score"]:
        if col not in markers.columns:
            markers[col] = np.nan

    markers["gene"] = markers["gene"].astype(str).str.strip()
    markers["cluster"] = markers["cluster"].astype(str).str.strip()
    numeric_cols = ["marker_logfc", "marker_p_adj", "marker_p_value", "marker_score"]
    for col in numeric_cols:
        markers[col] = pd.to_numeric(markers[col], errors="coerce")

    markers = markers.dropna(subset=["gene", "cluster", "marker_logfc", "marker_p_adj"])
    markers = markers[(markers["gene"] != "") & (markers["cluster"] != "")]
    markers = markers[
        (markers["marker_logfc"] > logfc_min)
        & (markers["marker_p_adj"] <= p_adj_cutoff)
    ]
    markers["marker_relevance"] = marker_relevance_score(
        markers["marker_logfc"],
        markers["marker_p_adj"],
    )

    markers = markers.sort_values(
        ["cluster", "gene", "marker_logfc", "marker_p_adj"],
        ascending=[True, True, False, True],
    )
    markers = markers.drop_duplicates(["cluster", "gene"], keep="first")

    grouped = (
        markers.groupby("cluster")["gene"]
        .apply(lambda values: sorted(set(values)))
        .to_dict()
    )
    if not grouped:
        raise ValueError(
            f"No full marker genes pass logFC > {logfc_min:g} and "
            f"adjusted p-value <= {p_adj_cutoff:g} in {marker_path}"
        )

    normalized = normalized_cluster_map(grouped, marker_path)
    marker_details = {}
    for cluster, cluster_markers in markers.groupby("cluster"):
        marker_details[cluster] = {
            row.gene: {
                "marker_logfc": row.marker_logfc,
                "marker_p_value": row.marker_p_value,
                "marker_p_adj": row.marker_p_adj,
                "marker_score": row.marker_score,
                "marker_relevance": row.marker_relevance,
            }
            for row in cluster_markers.itertuples(index=False)
        }

    return grouped, normalized, marker_details


def normalized_cluster_map(grouped, marker_path):
    normalized = {}
    for cluster in grouped:
        key = normalize_label(cluster)
        if key in normalized:
            raise ValueError(
                f"Marker clusters '{normalized[key]}' and '{cluster}' both normalize "
                f"to '{key}' in {marker_path}"
            )
        normalized[key] = cluster
    return normalized


def find_column(df, candidates, marker_path, required=True):
    for candidate in candidates:
        if candidate in df.columns:
            return candidate
    if required:
        raise ValueError(
            f"{marker_path} is missing one of the expected columns: {candidates}"
        )
    return None


def marker_relevance_score(logfc, p_adj):
    logfc = np.asarray(logfc, dtype=float)
    p_adj = np.asarray(p_adj, dtype=float)
    clipped_p = np.clip(p_adj, P_ADJ_FLOOR, 1.0)
    return np.maximum(logfc, 0.0) * -np.log10(clipped_p)


def discover_signature_paths(dataset_dir):
    paths = sorted(
        dataset_dir.glob("v*/signatures.csv"),
        key=lambda path: (parse_tau(path.parent.name), path.parent.name),
    )
    if not paths:
        raise FileNotFoundError(f"No v*/signatures.csv files found in {dataset_dir}")
    return paths


def parse_tau(version_name):
    if not version_name.startswith("v"):
        raise ValueError(f"Expected tau directory to start with 'v': {version_name}")
    return float(version_name[1:])


def format_tau(tau):
    return f"{tau:g}"


def benjamini_hochberg(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    qvalues = np.full(pvalues.shape, np.nan, dtype=float)
    valid = np.isfinite(pvalues)
    if not valid.any():
        return qvalues

    valid_indices = np.where(valid)[0]
    sorted_indices = valid_indices[np.argsort(pvalues[valid])]
    sorted_pvalues = pvalues[sorted_indices]
    m = len(sorted_pvalues)

    adjusted = sorted_pvalues * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    qvalues[sorted_indices] = np.minimum(adjusted, 1.0)
    return qvalues


def signature_column_map(signatures):
    normalized = {}
    for column in signatures.columns:
        key = normalize_label(column)
        if key in normalized:
            raise ValueError(
                f"Signature columns '{normalized[key]}' and '{column}' both normalize "
                f"to '{key}'"
            )
        normalized[key] = column
    return normalized


def analyze_signature_file(
    signature_path,
    marker_groups,
    marker_norm_to_cluster,
    marker_details,
    top_n,
):
    tau = parse_tau(signature_path.parent.name)
    signatures = pd.read_csv(signature_path, index_col=0)
    signatures.index = signatures.index.astype(str).str.strip()
    signatures = signatures.loc[signatures.index != ""]
    signatures = signatures.apply(pd.to_numeric, errors="coerce")

    universe = set(signatures.index)
    universe_size = len(universe)
    if universe_size == 0:
        raise ValueError(f"No genes found in {signature_path}")

    sig_norm_to_column = signature_column_map(signatures)

    missing_marker_clusters = sorted(set(marker_norm_to_cluster) - set(sig_norm_to_column))
    matched_marker_clusters = sorted(set(marker_norm_to_cluster) & set(sig_norm_to_column))

    summary_rows = []
    match_rows = []
    recovery_rows = []

    for norm_cluster in matched_marker_clusters:
        cluster = marker_norm_to_cluster[norm_cluster]
        column = sig_norm_to_column[norm_cluster]
        scores = signatures[column].dropna().sort_values(ascending=False)
        selected = scores.head(top_n)
        top_n_used = len(selected)

        markers = set(marker_groups[cluster])
        markers_in_universe = markers & universe
        top_genes = list(selected.index)
        overlapping_genes = set(top_genes) & markers_in_universe

        marker_count = len(markers_in_universe)
        overlap_count = len(overlapping_genes)
        precision = overlap_count / top_n_used if top_n_used else np.nan
        recall = overlap_count / marker_count if marker_count else np.nan

        if marker_count and top_n_used:
            pvalue = hypergeom.sf(
                overlap_count - 1,
                universe_size,
                marker_count,
                top_n_used,
            )
            expected_precision = marker_count / universe_size
            fold_enrichment = precision / expected_precision
        else:
            pvalue = np.nan
            fold_enrichment = np.nan

        continuous_metrics, cluster_recovery_rows = continuous_marker_metrics(
            scores=scores,
            marker_details_for_cluster=marker_details.get(cluster, {}),
            top_n_used=top_n_used,
            tau=tau,
            cell_type=cluster,
            signature_column=column,
        )
        recovery_rows.extend(cluster_recovery_rows)

        summary_rows.append(
            {
                "tau": tau,
                "tau_label": format_tau(tau),
                "cell_type": cluster,
                "signature_column": column,
                "universe_genes": universe_size,
                "top_n_requested": top_n,
                "top_n_used": top_n_used,
                "marker_genes_total": len(markers),
                "marker_genes_in_universe": marker_count,
                "overlap_count": overlap_count,
                "precision_at_n": precision,
                "recall": recall,
                "fold_enrichment": fold_enrichment,
                "p_value": pvalue,
                **continuous_metrics,
            }
        )

        for rank, (gene, score) in enumerate(selected.items(), start=1):
            gene_marker_details = marker_details.get(cluster, {}).get(gene, {})
            match_rows.append(
                {
                    "tau": tau,
                    "tau_label": format_tau(tau),
                    "cell_type": cluster,
                    "signature_column": column,
                    "rank": rank,
                    "gene": gene,
                    "signature_score": score,
                    "is_marker": gene in markers_in_universe,
                    "marker_logfc": gene_marker_details.get("marker_logfc", np.nan),
                    "marker_p_value": gene_marker_details.get("marker_p_value", np.nan),
                    "marker_p_adj": gene_marker_details.get("marker_p_adj", np.nan),
                    "marker_score": gene_marker_details.get("marker_score", np.nan),
                    "marker_relevance": gene_marker_details.get(
                        "marker_relevance",
                        np.nan,
                    ),
                }
            )

    match_info = {
        "matched": len(matched_marker_clusters),
        "total": len(marker_norm_to_cluster),
        "missing": [marker_norm_to_cluster[key] for key in missing_marker_clusters],
    }
    return summary_rows, match_rows, recovery_rows, match_info


def continuous_marker_metrics(
    scores,
    marker_details_for_cluster,
    top_n_used,
    tau,
    cell_type,
    signature_column,
):
    metric_names = {
        "marker_relevance_total_in_universe": np.nan,
        "marker_relevance_top_n": np.nan,
        "marker_relevance_ideal_top_n": np.nan,
        "weighted_precision_at_n": np.nan,
        "ndcg_at_n": np.nan,
        "recovery_auc_at_n": np.nan,
        "spearman_logfc": np.nan,
        "spearman_logfc_p_value": np.nan,
        "spearman_relevance": np.nan,
        "spearman_relevance_p_value": np.nan,
        "spearman_n": 0,
    }
    if not marker_details_for_cluster or top_n_used == 0:
        return metric_names, []

    marker_rows = []
    score_genes = set(scores.index)
    for gene, details in marker_details_for_cluster.items():
        if gene not in score_genes:
            continue
        marker_relevance = details.get("marker_relevance", np.nan)
        marker_logfc = details.get("marker_logfc", np.nan)
        if not np.isfinite(marker_relevance) or not np.isfinite(marker_logfc):
            continue
        marker_rows.append(
            {
                "gene": gene,
                "signature_score": scores.loc[gene],
                "marker_logfc": marker_logfc,
                "marker_relevance": marker_relevance,
            }
        )

    if not marker_rows:
        return metric_names, []

    marker_df = pd.DataFrame(marker_rows)
    marker_relevance_total = marker_df["marker_relevance"].sum()
    if not np.isfinite(marker_relevance_total) or marker_relevance_total <= 0:
        return metric_names, []

    top_genes = list(scores.index[:top_n_used])
    relevance_by_gene = dict(zip(marker_df["gene"], marker_df["marker_relevance"]))
    top_relevance = np.array(
        [relevance_by_gene.get(gene, 0.0) for gene in top_genes],
        dtype=float,
    )

    ideal_relevance = (
        marker_df["marker_relevance"]
        .sort_values(ascending=False)
        .head(top_n_used)
        .to_numpy(dtype=float)
    )
    if len(ideal_relevance) < top_n_used:
        ideal_relevance = np.pad(
            ideal_relevance,
            (0, top_n_used - len(ideal_relevance)),
            constant_values=0.0,
        )

    marker_relevance_top_n = top_relevance.sum()
    marker_relevance_ideal_top_n = ideal_relevance.sum()
    weighted_precision = (
        marker_relevance_top_n / marker_relevance_ideal_top_n
        if marker_relevance_ideal_top_n > 0
        else np.nan
    )
    ndcg = ndcg_score(top_relevance, ideal_relevance)
    recovery_auc = normalized_recovery_auc(
        top_relevance,
        ideal_relevance,
        marker_relevance_total,
    )

    spearman_logfc, spearman_logfc_p = safe_spearman(
        marker_df["signature_score"],
        marker_df["marker_logfc"],
    )
    spearman_relevance, spearman_relevance_p = safe_spearman(
        marker_df["signature_score"],
        marker_df["marker_relevance"],
    )

    metric_names.update(
        {
            "marker_relevance_total_in_universe": marker_relevance_total,
            "marker_relevance_top_n": marker_relevance_top_n,
            "marker_relevance_ideal_top_n": marker_relevance_ideal_top_n,
            "weighted_precision_at_n": weighted_precision,
            "ndcg_at_n": ndcg,
            "recovery_auc_at_n": recovery_auc,
            "spearman_logfc": spearman_logfc,
            "spearman_logfc_p_value": spearman_logfc_p,
            "spearman_relevance": spearman_relevance,
            "spearman_relevance_p_value": spearman_relevance_p,
            "spearman_n": len(marker_df),
        }
    )

    cumulative_relevance = np.cumsum(top_relevance)
    cumulative_ideal_relevance = np.cumsum(ideal_relevance)
    recovery_rows = []
    for rank in range(1, top_n_used + 1):
        recovery_rows.append(
            {
                "tau": tau,
                "tau_label": format_tau(tau),
                "cell_type": cell_type,
                "signature_column": signature_column,
                "rank": rank,
                "rank_fraction": rank / top_n_used,
                "recovered_relevance": cumulative_relevance[rank - 1],
                "ideal_relevance": cumulative_ideal_relevance[rank - 1],
                "recovered_relevance_fraction": (
                    cumulative_relevance[rank - 1] / marker_relevance_total
                ),
                "ideal_relevance_fraction": (
                    cumulative_ideal_relevance[rank - 1] / marker_relevance_total
                ),
            }
        )

    return metric_names, recovery_rows


def ndcg_score(top_relevance, ideal_relevance):
    discounts = 1 / np.log2(np.arange(2, len(top_relevance) + 2))
    dcg = np.sum(top_relevance * discounts)
    idcg = np.sum(ideal_relevance * discounts)
    return dcg / idcg if idcg > 0 else np.nan


def normalized_recovery_auc(top_relevance, ideal_relevance, total_relevance):
    if total_relevance <= 0 or len(top_relevance) == 0:
        return np.nan

    x = np.concatenate([[0.0], np.arange(1, len(top_relevance) + 1) / len(top_relevance)])
    recovered = np.concatenate([[0.0], np.cumsum(top_relevance) / total_relevance])
    ideal = np.concatenate([[0.0], np.cumsum(ideal_relevance) / total_relevance])
    ideal_auc = np.trapezoid(ideal, x)
    if ideal_auc <= 0:
        return np.nan
    return np.trapezoid(recovered, x) / ideal_auc


def safe_spearman(x, y):
    x = pd.Series(x, dtype=float)
    y = pd.Series(y, dtype=float)
    valid = x.notna() & y.notna()
    x = x[valid]
    y = y[valid]
    if len(x) < 3 or x.nunique() < 2 or y.nunique() < 2:
        return np.nan, np.nan
    rho, pvalue = spearmanr(x, y)
    return rho, pvalue


def analyze_dataset(
    dataset,
    results_path,
    marker_path,
    top_n,
    alpha,
    make_plots,
    marker_mode,
    output_name,
    marker_p_adj,
    marker_logfc_min,
):
    dataset_dir = results_path / dataset
    if marker_mode == "full":
        marker_groups, marker_norm_to_cluster, marker_details = read_full_marker_genes(
            marker_path,
            p_adj_cutoff=marker_p_adj,
            logfc_min=marker_logfc_min,
        )
    else:
        marker_groups, marker_norm_to_cluster, marker_details = read_plain_marker_genes(
            marker_path
        )
    signature_paths = discover_signature_paths(dataset_dir)

    all_summary_rows = []
    all_match_rows = []
    all_recovery_rows = []
    match_counts = []

    for signature_path in signature_paths:
        summary_rows, match_rows, recovery_rows, match_info = analyze_signature_file(
            signature_path,
            marker_groups,
            marker_norm_to_cluster,
            marker_details,
            top_n,
        )
        all_summary_rows.extend(summary_rows)
        all_match_rows.extend(match_rows)
        all_recovery_rows.extend(recovery_rows)
        match_counts.append(match_info)

    summary = pd.DataFrame(all_summary_rows)
    matches = pd.DataFrame(all_match_rows)
    recovery = pd.DataFrame(all_recovery_rows)
    if summary.empty:
        raise ValueError(f"No marker/signature matches could be analyzed for {dataset}")

    summary["q_value"] = benjamini_hochberg(summary["p_value"].values)
    summary["significant"] = summary["q_value"] < alpha
    summary["significance"] = summary["q_value"].apply(pvalue_to_stars)
    summary["marker_mode"] = marker_mode
    summary["marker_path"] = str(marker_path)
    summary["marker_p_adj_cutoff"] = marker_p_adj if marker_mode == "full" else np.nan
    summary["marker_logfc_min"] = marker_logfc_min if marker_mode == "full" else np.nan
    summary = summary.sort_values(["tau", "cell_type"]).reset_index(drop=True)
    matches["marker_mode"] = marker_mode
    matches["marker_path"] = str(marker_path)
    matches = matches.sort_values(["tau", "cell_type", "rank"]).reset_index(drop=True)
    if not recovery.empty:
        recovery["marker_mode"] = marker_mode
        recovery["marker_path"] = str(marker_path)
        recovery = recovery.sort_values(["tau", "cell_type", "rank"]).reset_index(
            drop=True
        )

    output_dir = dataset_dir / "signatures" / output_name
    output_dir.mkdir(parents=True, exist_ok=True)
    matches.to_csv(output_dir / "signature_marker_matches.csv", index=False)
    summary.to_csv(output_dir / "signature_marker_summary.csv", index=False)
    if not recovery.empty:
        recovery.to_csv(output_dir / "signature_marker_recovery_curves.csv", index=False)

    if marker_mode == "plain":
        matches.to_csv(dataset_dir / "signature_marker_matches.csv", index=False)
        summary.to_csv(dataset_dir / "signature_marker_summary.csv", index=False)

    if make_plots:
        plot_dataset(summary, recovery, output_dir, dataset, marker_mode)

    matched = min(item["matched"] for item in match_counts)
    total = max(item["total"] for item in match_counts)
    missing = sorted({name for item in match_counts for name in item["missing"]})

    return {
        "dataset": dataset,
        "marker_mode": marker_mode,
        "matched": matched,
        "total": total,
        "missing": missing,
        "n_tau": summary["tau"].nunique(),
        "n_tests": len(summary),
        "n_significant": int(summary["significant"].sum()),
    }


def plot_dataset(summary, recovery, plots_dir, dataset, marker_mode):
    plots_dir.mkdir(parents=True, exist_ok=True)

    sns.set_theme(style="whitegrid")
    sns.set_context(
        "talk",
        font_scale=1.1,
        rc={
            "axes.titlesize": 20,
            "axes.labelsize": 18,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize": 13,
        },
    )

    summary = summary.copy()
    tau_order = sorted(summary["tau"].unique())
    labels = [format_tau(tau) if tau < 1 else f"{format_tau(tau)} (NMF)" for tau in tau_order]
    marker_label = "Positive DE markers" if marker_mode == "full" else "Plain marker list"

    plot_metric_box(
        summary,
        tau_order,
        labels,
        "precision_at_n",
        "Precision@N",
        f"{dataset}: {marker_label} Recovery in Top Signature Genes",
        plots_dir / "signature_marker_precision",
    )

    plot_metric_box(
        summary,
        tau_order,
        labels,
        "fold_enrichment",
        "Fold Enrichment",
        f"{dataset}: {marker_label} Enrichment",
        plots_dir / "signature_marker_enrichment",
        reference_line=1.0,
    )

    plot_precision_heatmap(
        summary,
        tau_order,
        labels,
        f"{dataset}: {marker_label} Precision@N by Cell Type",
        plots_dir / "signature_marker_heatmap",
    )

    if marker_mode == "full":
        plot_full_marker_metrics(
            summary,
            recovery,
            tau_order,
            labels,
            dataset,
            marker_label,
            plots_dir,
        )


def plot_full_marker_metrics(
    summary,
    recovery,
    tau_order,
    labels,
    dataset,
    marker_label,
    plots_dir,
):
    plot_metric_box(
        summary,
        tau_order,
        labels,
        "ndcg_at_n",
        "NDCG@N",
        f"{dataset}: {marker_label} Ranking Quality",
        plots_dir / "signature_marker_ndcg",
        reference_line=1.0,
        annotate_enrichment=False,
    )

    plot_metric_box(
        summary,
        tau_order,
        labels,
        "weighted_precision_at_n",
        "Weighted Precision@N",
        f"{dataset}: {marker_label} Weighted Recovery",
        plots_dir / "signature_marker_weighted_precision",
        reference_line=1.0,
        annotate_enrichment=False,
    )

    plot_metric_box(
        summary,
        tau_order,
        labels,
        "recovery_auc_at_n",
        "Normalized Recovery AUC@N",
        f"{dataset}: {marker_label} Cumulative Recovery",
        plots_dir / "signature_marker_recovery_auc",
        reference_line=1.0,
        annotate_enrichment=False,
    )

    plot_metric_box(
        summary,
        tau_order,
        labels,
        "spearman_logfc",
        "Spearman rho",
        f"{dataset}: Signature Score vs Marker logFC",
        plots_dir / "signature_marker_spearman_logfc",
        reference_line=0.0,
        annotate_enrichment=False,
    )

    plot_metric_box(
        summary,
        tau_order,
        labels,
        "spearman_relevance",
        "Spearman rho",
        f"{dataset}: Signature Score vs Marker Relevance",
        plots_dir / "signature_marker_spearman_relevance",
        reference_line=0.0,
        annotate_enrichment=False,
    )

    plot_metric_heatmap(
        summary,
        tau_order,
        labels,
        "ndcg_at_n",
        "NDCG@N",
        f"{dataset}: {marker_label} NDCG@N by Cell Type",
        plots_dir / "signature_marker_ndcg_heatmap",
    )

    plot_metric_heatmap(
        summary,
        tau_order,
        labels,
        "spearman_logfc",
        "Spearman rho",
        f"{dataset}: Signature Score vs logFC by Cell Type",
        plots_dir / "signature_marker_spearman_logfc_heatmap",
        center=0.0,
        cmap="vlag",
    )

    if not recovery.empty:
        plot_recovery_curves(
            recovery,
            tau_order,
            labels,
            f"{dataset}: Mean Cumulative Marker Relevance Recovery",
            plots_dir / "signature_marker_recovery_curves",
        )


def plot_metric_box(
    summary,
    tau_order,
    labels,
    metric,
    ylabel,
    title,
    output_stem,
    reference_line=None,
    annotate_enrichment=True,
):
    plot_data = summary[np.isfinite(summary[metric])].copy()
    if plot_data.empty:
        return

    stats = (
        plot_data.groupby("tau", as_index=False)[metric]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    stats = stats.set_index("tau").reindex(tau_order).reset_index()
    stats["sem"] = stats["std"] / np.sqrt(stats["count"])
    stats["sem"] = stats["sem"].fillna(0.0)

    x_positions = np.arange(len(tau_order))
    x_lookup = dict(zip(tau_order, x_positions))
    plot_data["x_position"] = plot_data["tau"].map(x_lookup)

    plt.figure(figsize=(11, 6))
    color = sns.color_palette("viridis", n_colors=len(tau_order))[max(0, len(tau_order) // 2)]
    ax = plt.gca()
    sns.stripplot(
        data=plot_data,
        x="x_position",
        y=metric,
        order=x_positions,
        color="black",
        alpha=0.28,
        size=3.5,
        jitter=0.16,
        ax=ax,
    )
    ax.errorbar(
        x_positions,
        stats["mean"],
        yerr=stats["sem"],
        fmt="o-",
        color=color,
        ecolor=color,
        elinewidth=2,
        capsize=5,
        capthick=2,
        markersize=7,
        linewidth=2.4,
        zorder=10,
    )
    if reference_line is not None:
        ax.axhline(reference_line, color="0.3", linewidth=1.2, linestyle="--")
    if annotate_enrichment:
        add_tau_significance_counts(ax, plot_data, tau_order, metric)
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel(ylabel)
    ax.set_xticks(range(len(tau_order)))
    ax.set_xticklabels(labels, rotation=30, ha="right")
    plt.tight_layout()
    save_plot(output_stem)


def pvalue_to_stars(pvalue):
    if not np.isfinite(pvalue):
        return ""
    if pvalue < 1e-4:
        return "****"
    if pvalue < 1e-3:
        return "***"
    if pvalue < 1e-2:
        return "**"
    if pvalue < 0.05:
        return "*"
    return ""


def add_tau_significance_counts(ax, summary, tau_order, metric):
    ymin, ymax = ax.get_ylim()
    y_range = ymax - ymin
    if not np.isfinite(y_range) or y_range <= 0:
        y_range = max(abs(summary[metric].max()), 1.0)

    gap = y_range * 0.06
    top = ymax

    for index, tau in enumerate(tau_order):
        values = summary[np.isclose(summary["tau"], tau)]
        significant = values["significance"].astype(bool).sum()
        total = len(values)
        if significant == 0 or total == 0:
            label = "ns"
        else:
            strongest = min(values["q_value"].dropna())
            label = f"{pvalue_to_stars(strongest)}\n{significant}/{total}"

        y = values[metric].max() + gap
        ax.text(
            index,
            y,
            label,
            ha="center",
            va="bottom",
            fontsize=11,
            linespacing=0.85,
        )
        top = max(top, y + gap)

    ax.set_ylim(ymin, top)


def plot_precision_heatmap(summary, tau_order, labels, title, output_stem):
    heatmap_data = summary.pivot(index="cell_type", columns="tau", values="precision_at_n")
    heatmap_data = heatmap_data.reindex(columns=tau_order)
    heatmap_data = heatmap_data.loc[heatmap_data.mean(axis=1).sort_values(ascending=False).index]
    star_data = summary.pivot(index="cell_type", columns="tau", values="significance")
    star_data = star_data.reindex(index=heatmap_data.index, columns=tau_order).fillna("")
    annotations = heatmap_data.map(lambda value: f"{value:.2f}" if np.isfinite(value) else "")
    annotations = annotations + star_data.map(lambda value: f"\n{value}" if value else "")

    height = max(6, min(16, 0.35 * len(heatmap_data) + 2))
    plt.figure(figsize=(12, height))
    ax = sns.heatmap(
        heatmap_data,
        annot=annotations,
        fmt="",
        cmap="viridis",
        vmin=0,
        vmax=max(float(np.nanmax(heatmap_data.values)), 1e-8),
        linewidths=0.3,
        linecolor="white",
        cbar_kws={"label": "Precision@N"},
    )
    ax.set_title(title)
    ax.set_xlabel("tau")
    ax.set_ylabel("")
    ax.set_xticks(np.arange(len(tau_order)) + 0.5)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    plt.tight_layout()
    save_plot(output_stem)


def plot_metric_heatmap(
    summary,
    tau_order,
    labels,
    metric,
    cbar_label,
    title,
    output_stem,
    center=None,
    cmap="viridis",
):
    heatmap_data = summary.pivot(index="cell_type", columns="tau", values=metric)
    heatmap_data = heatmap_data.reindex(columns=tau_order)
    heatmap_data = heatmap_data.dropna(axis=0, how="all")
    if heatmap_data.empty:
        return

    heatmap_data = heatmap_data.loc[
        heatmap_data.mean(axis=1).sort_values(ascending=False).index
    ]
    annotations = heatmap_data.map(
        lambda value: f"{value:.2f}" if np.isfinite(value) else ""
    )

    height = max(6, min(16, 0.35 * len(heatmap_data) + 2))
    plt.figure(figsize=(12, height))
    heatmap_kwargs = {
        "data": heatmap_data,
        "annot": annotations,
        "fmt": "",
        "cmap": cmap,
        "linewidths": 0.3,
        "linecolor": "white",
        "cbar_kws": {"label": cbar_label},
    }
    if center is None:
        heatmap_kwargs["vmin"] = 0
        heatmap_kwargs["vmax"] = max(float(np.nanmax(heatmap_data.values)), 1e-8)
    else:
        heatmap_kwargs["center"] = center

    ax = sns.heatmap(**heatmap_kwargs)
    ax.set_title(title)
    ax.set_xlabel("tau")
    ax.set_ylabel("")
    ax.set_xticks(np.arange(len(tau_order)) + 0.5)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    plt.tight_layout()
    save_plot(output_stem)


def plot_recovery_curves(recovery, tau_order, labels, title, output_stem):
    plot_data = (
        recovery.groupby(["tau", "rank"], as_index=False)
        .agg(
            recovered_relevance_fraction=("recovered_relevance_fraction", "mean"),
            ideal_relevance_fraction=("ideal_relevance_fraction", "mean"),
        )
        .sort_values(["tau", "rank"])
    )
    if plot_data.empty:
        return

    label_map = dict(zip(tau_order, labels))
    plot_data["tau_label_plot"] = plot_data["tau"].map(label_map)

    plt.figure(figsize=(11, 7))
    ax = sns.lineplot(
        data=plot_data,
        x="rank",
        y="recovered_relevance_fraction",
        hue="tau_label_plot",
        hue_order=labels,
        palette=sns.color_palette("viridis", n_colors=len(tau_order)),
        linewidth=2,
    )
    ax.set_title(title)
    ax.set_xlabel("Top signature genes")
    ax.set_ylabel("Recovered marker relevance fraction")
    ax.legend(title="tau", ncols=2, fontsize=10, title_fontsize=11)
    plt.tight_layout()
    save_plot(output_stem)


def save_plot(output_stem):
    plt.savefig(f"{output_stem}.png", dpi=300)
    plt.savefig(f"{output_stem}.pdf", dpi=300)
    plt.close()


def main():
    args = parse_args()
    if args.top_n <= 0:
        raise ValueError("--top-n must be positive")
    if not 0 < args.alpha < 1:
        raise ValueError("--alpha must be between 0 and 1")
    if not 0 <= args.marker_p_adj <= 1:
        raise ValueError("--marker-p-adj must be between 0 and 1")

    marker_modes = (
        ["plain", "full"] if args.marker_mode == "both" else [args.marker_mode]
    )
    results = []
    for marker_mode in marker_modes:
        marker_paths = (
            DEFAULT_FULL_MARKER_PATHS if marker_mode == "full" else DEFAULT_MARKER_PATHS
        )
        output_name = "full_dataframe" if marker_mode == "full" else "plain_list"

        if args.dataset == "all":
            datasets = marker_paths.keys()
        elif args.dataset in marker_paths:
            datasets = [args.dataset]
        elif args.marker_mode == "both" and marker_mode == "full":
            continue
        else:
            available = ", ".join(marker_paths)
            raise ValueError(
                f"{marker_mode!r} marker mode is not configured for {args.dataset}. "
                f"Available datasets: {available}"
            )

        for dataset in datasets:
            result = analyze_dataset(
                dataset=dataset,
                results_path=args.results_path,
                marker_path=marker_paths[dataset],
                top_n=args.top_n,
                alpha=args.alpha,
                make_plots=not args.no_plots,
                marker_mode=marker_mode,
                output_name=output_name,
                marker_p_adj=args.marker_p_adj,
                marker_logfc_min=args.marker_logfc_min,
            )
            results.append(result)

    for result in results:
        status = f"{result['matched']}/{result['total']}"
        line = (
            f"{result['dataset']} ({result['marker_mode']}): "
            f"matched marker clusters {status}; "
            f"tau values={result['n_tau']}; tests={result['n_tests']}; "
            f"significant(q<{args.alpha:g})={result['n_significant']}"
        )
        if result["missing"]:
            line += "; missing=" + ", ".join(result["missing"])
        print(line)


if __name__ == "__main__":
    main()
