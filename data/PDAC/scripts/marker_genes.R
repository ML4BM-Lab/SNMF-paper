
script_arg <- grep("--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("--file=", "", script_arg[1])))
} else {
  getwd()
}
data_dir <- normalizePath(file.path(script_dir, ".."))
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."))
zenodo_dir <- file.path(repo_root, "data", "zenodo", "PDAC")
dir.create(zenodo_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(data_dir, "sc_count.RData"))
load(file.path(data_dir, "sc_meta.RData"))

library(dplyr)
library(Seurat)

# Create Seurat object
sc <- CreateSeuratObject(counts = sc_count, meta.data = sc_meta)

Idents(sc) <- sc$cellType

# Normalization & scaling (needed for marker detection)
sc <- NormalizeData(sc)
sc <- FindVariableFeatures(sc)
sc <- ScaleData(sc)

all_markers <- FindAllMarkers(
  sc, 
  only.pos = FALSE,     
  min.pct = 0.25,          # expressed in ≥25% of cells
  logfc.threshold = 0.25   # minimum log2 fold-change
)

sig_markers <- all_markers %>% filter(p_val_adj < 0.05)
write.csv(sig_markers, file.path(zenodo_dir, "PDAC_marker_genes_full.csv"), row.names = FALSE)

# convert to list by cluster
marker_list <- split(sig_markers$gene, sig_markers$cluster)

marker_df <- stack(marker_list)
colnames(marker_df) <- c("gene", "cluster")

write.csv(marker_df, file.path(zenodo_dir, "PDAC_marker_genes.csv"), row.names = FALSE)
