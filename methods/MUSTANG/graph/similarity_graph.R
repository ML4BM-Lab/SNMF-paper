#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
output_path <- args[1]

transcrip_neighbors <- read.csv(paste0(output_path, "transcrip_adjacency.csv"))
spatial_neighbors   <- read.csv(paste0(output_path, "spatial_neighbors.csv"))
colnames(spatial_neighbors) <- colnames(transcrip_neighbors)

similarity_neighbors <- rbind(transcrip_neighbors, spatial_neighbors)
write.csv(similarity_neighbors, paste0(output_path, "similarity_neighbors.csv"), quote=FALSE, row.names=TRUE)
