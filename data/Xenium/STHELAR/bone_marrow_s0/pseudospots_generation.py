
from pathlib import Path

import numpy as np
import pandas as pd

import scanpy as sc
import anndata as ad
from scipy import sparse

from sklearn.neighbors import KDTree

# =========================================================
# USER CONFIGURATION
# =========================================================

SDATA_PATH = Path(
    "/scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/bone_marrow_s0/sdata_bone_marrow_s0.zarr"
)

TABLE_NAME = "table_cells"
CELL_TYPE_LABEL = "final_label"
QCFILTER_CELLTYPE_LABEL = "Less10"
COUNTS_LAYER_NAME = "count"

SPOT_DIAMETER_UM = 55.0
SPOT_RADIUS_UM = SPOT_DIAMETER_UM / 2.0
CENTER_TO_CENTER_UM = 100.0

MIN_CELLS_PER_SPOT = 1

OUT_DIR = Path("/scratch/lalonsoeste/PhD/NMF_deconvolution/data/Xenium/STHELAR/bone_marrow_s0") / "pseudospots"
OUT_H5AD = OUT_DIR / "pseudospots.h5ad"
OUT_PROPS = OUT_DIR / "proportions.csv"
OUT_MARKER_GENES = OUT_DIR / "marker_genes.csv"

# =========================================================
# LOAD TABLE_CELLS DIRECTLY
# =========================================================

TABLE_PATH = SDATA_PATH / "tables" / TABLE_NAME
adata = ad.read_zarr(TABLE_PATH)

print("\nLoaded AnnData:")
print(adata)

print("\nobs columns:")
print(list(adata.obs.columns))

print("\nobsm keys:")
print(list(adata.obsm.keys()))

if CELL_TYPE_LABEL not in adata.obs.columns:
    raise ValueError(f"{CELL_TYPE_LABEL!r} not found in adata.obs")
else:
    adata = adata[adata.obs[CELL_TYPE_LABEL] != QCFILTER_CELLTYPE_LABEL,:].copy()

# =========================================================
# GET COORDINATES
# =========================================================

if "spatial" in adata.obsm:
    coords = np.asarray(adata.obsm["spatial"])
    print("Using coordinates from adata.obsm['spatial']")
else:
    coord_candidates = [
        ("x", "y"),
        ("center_x", "center_y"),
        ("centroid_x", "centroid_y"),
        ("cell_centroid_x", "cell_centroid_y"),
        ("global_x", "global_y"),
    ]

    coords = None
    for xcol, ycol in coord_candidates:
        if xcol in adata.obs.columns and ycol in adata.obs.columns:
            coords = adata.obs[[xcol, ycol]].to_numpy()
            print(f"Using coordinates from adata.obs: {xcol}, {ycol}")
            break

    if coords is None:
        raise ValueError("Could not find spatial coordinates.")

coords = coords[:, :2].astype(float)

# =========================================================
# LABELS + EXPRESSION
# =========================================================

labels = adata.obs[CELL_TYPE_LABEL].astype(str).to_numpy()
cell_types = np.sort(pd.unique(labels))

X = adata.layers[COUNTS_LAYER_NAME]

# =========================================================
# BUILD VISIUM-LIKE HEXAGONAL GRID
# =========================================================

x = coords[:, 0]
y = coords[:, 1]

dx = CENTER_TO_CENTER_UM
dy = CENTER_TO_CENTER_UM * np.sqrt(3) / 2.0

xmin, xmax = x.min(), x.max()
ymin, ymax = y.min(), y.max()

centers = []
grid_coords = []

row = 0
yy = ymin

while yy <= ymax:
    x_offset = 0.0 if row % 2 == 0 else dx / 2.0
    xx = xmin + x_offset
    col = int(x_offset > 0.0)

    while xx <= xmax:
        centers.append((xx, yy))
        xx += dx

        grid_coords.append((col, row))
        col += 2

    yy += dy
    row += 1

centers = np.asarray(centers)
grid_coords = np.asarray(grid_coords)

print(f"\nCandidate Visium-like spots: {len(centers)}")

# =========================================================
# FIND CELLS INSIDE EACH 55 um SPOT
# =========================================================

tree = KDTree(coords)

spot_to_cells = tree.query_radius(
    centers,
    r=SPOT_RADIUS_UM
)

valid = np.array([
    len(idx) >= MIN_CELLS_PER_SPOT
    for idx in spot_to_cells
])

centers = centers[valid]
grid_coords = grid_coords[valid]
spot_to_cells = spot_to_cells[valid]

print(f"Non-empty spots: {len(centers)}")

# =========================================================
# POOL COUNTS + GROUND-TRUTH PROPORTIONS
# =========================================================

spot_expr = []
spot_props = []
spot_ncells = []

for idx in spot_to_cells:
    Xi = X[idx]

    if sparse.issparse(Xi):
        expr_sum = np.asarray(Xi.sum(axis=0)).ravel()
    else:
        expr_sum = np.asarray(Xi.sum(axis=0)).ravel()

    spot_expr.append(expr_sum)

    spot_labels = labels[idx]
    counts = pd.Series(spot_labels).value_counts()

    props = {
        ct: counts.get(ct, 0) / len(idx)
        for ct in cell_types
    }

    spot_props.append(props)
    spot_ncells.append(len(idx))

spot_expr = np.vstack(spot_expr)
spot_props = pd.DataFrame(spot_props)
spot_props = spot_props.loc[:,spot_props.sum(axis=0) > 0]
print(f"Number of non-zero cell types: {spot_props.shape[1]}")

# =========================================================
# MARKER GENES
# =========================================================

sc.tl.rank_genes_groups(
    adata,
    layer="log_norm",
    groupby=CELL_TYPE_LABEL,
    method="wilcoxon"
)

result = adata.uns["rank_genes_groups"]
groups = result["names"].dtype.names

dfs = []

for group in groups:
    df = pd.DataFrame({
        "gene": result["names"][group],
        "logfoldchange": result["logfoldchanges"][group],
        "pvals_adj": result["pvals_adj"][group],
        "score": result["scores"][group]
    })

    df = df[
        (df["pvals_adj"] < 0.05) &
        (df["logfoldchange"] > 0.25)
    ]

    df["cluster"] = group
    dfs.append(df)

markers_df = pd.concat(dfs, ignore_index=True)

markers_df = markers_df[
    ["gene", "cluster"]
]


# =========================================================
# SAVE OUTPUTS
# =========================================================

OUT_DIR.mkdir(parents=True, exist_ok=True)

spot_ids = [
    f"{x}x{y}"
    for x,y in grid_coords
]

adata_spots = ad.AnnData(
    X=sparse.csr_matrix(spot_expr),
    obs=pd.DataFrame(index=spot_ids),
    var=adata.var.copy(),
)

adata_spots.obsm["spatial"] = centers

adata_spots.obs["n_cells"] = spot_ncells
adata_spots.obs["spot_diameter_um"] = SPOT_DIAMETER_UM
adata_spots.obs["spot_radius_um"] = SPOT_RADIUS_UM
adata_spots.obs["center_to_center_um"] = CENTER_TO_CENTER_UM
adata_spots.obs["source_sdata"] = str(SDATA_PATH)
adata_spots.obs["source_table"] = TABLE_NAME
adata_spots.obs["label_column_used"] = CELL_TYPE_LABEL

sc.pp.filter_cells(adata_spots, min_counts=1)
sc.pp.filter_genes(adata_spots, min_counts=1)

spot_props.index = spot_ids
spot_props = spot_props.loc[adata_spots.obs_names,:]

adata_spots.obs = pd.concat([adata_spots.obs, spot_props], axis=1)

adata_spots.write_h5ad(OUT_H5AD)
adata_spots.to_df().transpose().to_csv(OUT_DIR / "pseudospots.csv")
spot_props.to_csv(OUT_PROPS)
markers_df.to_csv(OUT_MARKER_GENES)

sc.pl.spatial(adata_spots, color=spot_props.columns, spot_size=SPOT_DIAMETER_UM, show=False)[0].figure.savefig(OUT_DIR / "spatial_proportions.png")

print(f"\nSaved h5ad: {OUT_H5AD}")
print(f"Saved proportions: {OUT_PROPS}")
print(f"Saved marker genes: {OUT_MARKER_GENES}")

print("\nOutput AnnData:")
print(adata_spots)

print("\nProportions preview:")
print(spot_props.head())

print("\nMarker genes preview:")
print(markers_df.head())