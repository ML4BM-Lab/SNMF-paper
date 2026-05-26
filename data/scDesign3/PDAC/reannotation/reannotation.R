
library(dplyr)
library(jsonlite)

proportions <- read.csv("/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/counts/proportions.csv", row.names=1, check.names=FALSE)
rowSums(proportions)

marker_genes <- read.csv("/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/marker_genes/sig_marker_genes.csv", check.names=FALSE)

reannotate_lvl5 <- c(
  "Acinar_cells" = "Epithelial/Cancer",
  "Ductal_terminal_ductal_like" = "Epithelial/Cancer",
  "Ductal_CRISP3_high-centroacinar_like" = "Epithelial/Cancer",
  "Ductal_MHC_Class_II" = "Epithelial/Cancer",
  "Ductal_APOL1_high-hypoxic" = "Epithelial/Cancer",
  "Cancer_clone_A" = "Epithelial/Cancer",
  "Cancer_clone_B" = "Epithelial/Cancer",
  "Tuft_cells" = "Epithelial/Cancer",
  "Endocrine_cells" = "Endocrine",
  "mDCs_A" = "Immune",
  "mDCs_B" = "Immune",
  "pDCs" = "Immune",
  "Macrophages_A" = "Immune",
  "Macrophages_B" = "Immune",
  "Monocytes" = "Immune",
  "T_cells_&_NK_cells" = "Immune",
  "Mast_cells" = "Immune",
  "Fibroblasts" = "Stromal",
  "Endothelial_cells" = "Stromal",
  "RBCs" = "Other"
)

proportions_lvl5 <- as.data.frame(
  sapply(split(colnames(proportions), reannotate_lvl5[colnames(proportions)]),
         function(cols) rowSums(proportions[, cols, drop = FALSE]))
)
rowSums(proportions_lvl5)
write.csv(proportions_lvl5, file="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/reannotation/k5/proportions_k5.csv")

markers_lvl5 <- marker_genes %>%
  mutate(cluster = reannotate_lvl5[cluster])

write.csv(markers_lvl5, "/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/reannotation/k5/marker_genes_k5.csv")
write_json(reannotate_lvl5, "/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/reannotation/k5/mapping_k5.json", pretty = TRUE, auto_unbox = TRUE)

reannotate_lvl10 <- c(
  "Acinar_cells" = "Acinar",
  "Ductal_terminal_ductal_like" = "Ductal",
  "Ductal_CRISP3_high-centroacinar_like" = "Ductal",
  "Ductal_MHC_Class_II" = "Ductal",
  "Ductal_APOL1_high-hypoxic" = "Ductal",
  "Cancer_clone_A" = "Cancer",
  "Cancer_clone_B" = "Cancer",
  "Tuft_cells" = "Tuft",
  "Endocrine_cells" = "Endocrine",
  "Macrophages_A" = "Myeloid",
  "Macrophages_B" = "Myeloid",
  "Monocytes" = "Myeloid",
  "mDCs_A" = "Dendritic",
  "mDCs_B" = "Dendritic",
  "pDCs" = "Dendritic",
  "T_cells_&_NK_cells" = "Lymphoid",
  "Mast_cells" = "Mast",
  "Fibroblasts" = "Stromal",
  "Endothelial_cells" = "Stromal",
  "RBCs" = "Stromal"
)

proportions_lvl10 <- as.data.frame(
  sapply(split(colnames(proportions), reannotate_lvl10[colnames(proportions)]),
         function(cols) rowSums(proportions[, cols, drop = FALSE]))
)
rowSums(proportions_lvl10)
write.csv(proportions_lvl10, file="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/reannotation/proportions_k10.csv")

markers_lvl10 <- marker_genes %>%
  mutate(cluster = reannotate_lvl10[cluster])

write.csv(markers_lvl10, "/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/reannotation/marker_genes_k10.csv")
write_json(reannotate_lvl10, "/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/reannotation/mapping_k10.json", pretty = TRUE, auto_unbox = TRUE)

reannotate_lvl15 <- c(
  "Acinar_cells" = "Acinar",
  "Ductal_terminal_ductal_like" = "Ductal_terminal",
  "Ductal_CRISP3_high-centroacinar_like" = "Ductal_CRISP3_high",
  "Ductal_MHC_Class_II" = "Ductal_MHC_Class_II",
  "Ductal_APOL1_high-hypoxic" = "Ductal_APOL1_high",
  "Cancer_clone_A" = "Cancer_clone_A",
  "Cancer_clone_B" = "Cancer_clone_B",
  "Tuft_cells" = "Tuft",
  "Endocrine_cells" = "Endocrine",
  "Macrophages_A" = "Macrophage/Monocyte",
  "Macrophages_B" = "Macrophage/Monocyte",
  "Monocytes" = "Macrophage/Monocyte",
  "mDCs_A" = "mDCs",
  "mDCs_B" = "mDCs",
  "pDCs" = "pDCs",
  "T_cells_&_NK_cells" = "T/NK",
  "Mast_cells" = "Mast",
  "Fibroblasts" = "Stromal",
  "Endothelial_cells" = "Stromal",
  "RBCs" = "Stromal"
)

proportions_lvl15 <- as.data.frame(
  sapply(split(colnames(proportions), reannotate_lvl15[colnames(proportions)]),
         function(cols) rowSums(proportions[, cols, drop = FALSE]))
)
rowSums(proportions_lvl15)
write.csv(proportions_lvl15, file="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/reannotation/proportions_k15.csv")

markers_lvl15 <- marker_genes %>%
  mutate(cluster = reannotate_lvl15[cluster])

write.csv(markers_lvl15, "/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/reannotation/marker_genes_k15.csv")
write_json(reannotate_lvl15, "/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/reannotation/mapping_k15.json", pretty = TRUE, auto_unbox = TRUE)

reannotate_lvl8 <- c(
  "Acinar_cells" = "Acinar",
  "Ductal_terminal_ductal_like" = "Ductal",
  "Ductal_CRISP3_high-centroacinar_like" = "Ductal",
  "Ductal_MHC_Class_II" = "Ductal",
  "Ductal_APOL1_high-hypoxic" = "Ductal",
  "Cancer_clone_A" = "Cancer",
  "Cancer_clone_B" = "Cancer",
  "Tuft_cells" = "Tuft",
  "Endocrine_cells" = "Endocrine",
  "mDCs_A" = "Immune",
  "mDCs_B" = "Immune",
  "pDCs" = "Immune",
  "Macrophages_A" = "Immune",
  "Macrophages_B" = "Immune",
  "Monocytes" = "Immune",
  "T_cells_&_NK_cells" = "Immune",
  "Mast_cells" = "Immune",
  "Fibroblasts" = "Stromal",
  "Endothelial_cells" = "Stromal",
  "RBCs" = "RBCs"
)

proportions_lvl8 <- as.data.frame(
  sapply(split(colnames(proportions), reannotate_lvl8[colnames(proportions)]),
         function(cols) rowSums(proportions[, cols, drop = FALSE]))
)
rowSums(proportions_lvl8)
write.csv(proportions_lvl8, file="/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/reannotation/proportions_k8.csv")

markers_lvl8 <- marker_genes %>%
  mutate(cluster = reannotate_lvl8[cluster])

write.csv(markers_lvl8, "/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/reannotation/marker_genes_k8.csv")
write_json(reannotate_lvl8, "/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/reannotation/mapping_k8.json", pretty = TRUE, auto_unbox = TRUE)