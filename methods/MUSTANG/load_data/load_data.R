#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
data_path <- args[1]
output_path <- args[2]

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(BayesSpace)
  library(scater)
})

set.seed(100)

counts <- read.csv(data_path, row.names = 1, check.names = FALSE)

# Convert to matrix (or sparse matrix if large)
counts <- as.matrix(counts)

# --- If you have spatial coordinates ---
# For example, from the column names "10x13", split into x and y
locations <- data.frame(
  spot = colnames(counts)
)
coords <- do.call(rbind, strsplit(locations$spot, "x"))
coords <- data.frame(x = as.numeric(coords[,1]), y = as.numeric(coords[,2]))
rownames(coords) <- colnames(counts)

# --- Create SingleCellExperiment ---
sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = coords
)

sce <- spatialPreprocess(sce, n.PCs=50, n.HVGs=2000, assay.type="logcounts")
sce <- runUMAP(sce, dimred="PCA")
colnames(reducedDim(sce, "UMAP")) <- c("UMAP1", "UMAP2")

saveRDS(sce, file = paste0(output_path, "sce.rds"))
