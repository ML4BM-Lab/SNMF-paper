
load("./../sc_count.RData")
load("./../sc_meta.RData")

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

# convert to list by cluster
marker_list <- split(sig_markers$gene, sig_markers$cluster)

marker_df <- stack(marker_list)
colnames(marker_df) <- c("gene", "cluster")

write.csv(marker_df, "./marker_genes.csv", row.names = FALSE)


