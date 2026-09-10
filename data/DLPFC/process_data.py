#!/usr/bin/env python3

import argparse
import shutil
from pathlib import Path
from typing import List

import pandas as pd
import scanpy as sc
import tqdm


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "final"


def sample_id_from_dir(sample_dir: Path) -> str:
    return sample_dir.name


def existing_path(*candidates: Path) -> Path:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    expected = ", ".join(str(path) for path in candidates)
    raise FileNotFoundError(f"None of these input files exist: {expected}")


def process_sample(sample_dir: Path, output_root: Path) -> None:
    sample = sample_id_from_dir(sample_dir)
    output_dir = output_root / sample
    output_dir.mkdir(parents=True, exist_ok=True)

    h5ad_path = existing_path(
        sample_dir / "filtered_feature_bc_matrix.h5ad",
        sample_dir / f"{sample}_filtered_feature_bc_matrix.h5ad",
    )
    truth_path = existing_path(
        sample_dir / f"{sample}_truth.txt",
        sample_dir / "truth.txt",
    )

    final_h5ad = output_dir / f"{sample}_filtered_feature_bc_matrix.h5ad"
    final_truth = output_dir / f"{sample}_truth.txt"
    final_counts = output_dir / f"{sample}_counts.csv"
    final_markers = output_dir / f"{sample}_marker_genes.csv"
    final_k = output_dir / f"{sample}_K.txt"

    if h5ad_path.resolve() != final_h5ad.resolve():
        shutil.copy2(h5ad_path, final_h5ad)
    if truth_path.resolve() != final_truth.resolve():
        shutil.copy2(truth_path, final_truth)

    adata = sc.read_h5ad(final_h5ad)

    n_genes = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=int(adata.n_obs * 0.1))
    print(f"{sample}: filtered out {n_genes - adata.n_vars} genes; kept {adata.n_vars}")

    counts = adata.X.todense() if hasattr(adata.X, "todense") else adata.X
    pd.DataFrame(
        counts.T,
        index=adata.var_names,
        columns=[
            f"{spot.obs['array_row'].values[0]}x{spot.obs['array_col'].values[0]}"
            for spot in adata
        ],
    ).to_csv(final_counts)

    if "Region" not in adata.obs.columns:
        print(f"{sample}: assigning ground-truth regions")
        with final_truth.open() as handle:
            for line in tqdm.tqdm(handle):
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 2:
                    continue
                spot, cluster = parts
                if spot in adata.obs.index:
                    adata.obs.loc[spot, "Region"] = cluster

    k = adata.obs["Region"].nunique()
    final_k.write_text(f"{k}\n")
    print(f"{sample}: detected k={k}")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.tl.rank_genes_groups(adata, groupby="Region", method="wilcoxon")

    markers = sc.get.rank_genes_groups_df(adata, None)
    threshold = 0.05
    result = pd.DataFrame(columns=["cluster", "gene"])
    clusters_with_markers = -1

    while clusters_with_markers < k:
        if threshold >= 1:
            print(f"{sample}: WARNING, not all clusters are present in marker genes")
            break

        sig = markers[
            (markers["pvals_adj"] < threshold)
            & (markers["logfoldchanges"] > 0)
        ]
        result = sig[["group", "names"]].copy()
        result.columns = ["cluster", "gene"]
        clusters_with_markers = len(pd.unique(result["cluster"]))
        threshold += 0.01

    print(f"{sample}: marker genes cover {clusters_with_markers} clusters")
    result.to_csv(final_markers, index=False)


def iter_sample_dirs(source_dir: Path) -> List[Path]:
    sample_dirs = [path for path in sorted(source_dir.iterdir()) if path.is_dir()]
    if sample_dirs:
        return sample_dirs
    return [source_dir]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate final Zenodo-style DLPFC files from raw per-sample folders."
    )
    parser.add_argument(
        "source",
        type=Path,
        help="Raw DLPFC root containing sample folders, or a single raw sample folder.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Output root for final sample folders. Defaults to data/DLPFC/final.",
    )
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    for sample_dir in iter_sample_dirs(args.source.resolve()):
        process_sample(sample_dir, output_dir)


if __name__ == "__main__":
    main()
