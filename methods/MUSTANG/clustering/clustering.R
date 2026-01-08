#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
output_path <- args[1]

suppressPackageStartupMessages({
  library(MERINGUE)
  library(igraph)
})

sce <- readRDS(paste0(output_path, "sce.rds"))

k <- 20
clusters <- MERINGUE::getClusters(t(logcounts(sce)), k, weight=TRUE, method=igraph::cluster_louvain)
colData(sce)$Cluster <- clusters

saveRDS(sce, paste0(output_path, "clusters.rds"))
