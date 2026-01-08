#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
output_path <- args[1]

suppressPackageStartupMessages({library(SpatialExperiment)})

sce <- readRDS(paste0(output_path, "clusters.rds"))
meta <- data.frame(Cluster=colData(sce)$Cluster)

transcrip_edge_num <- 0
transcrip_adjacency <- list(Node1=integer(), Node2=integer())

for(i in 1:(nrow(meta)-1)){
  for(j in (i+1):nrow(meta)){
    if((meta$Cluster[i]==meta$Cluster[j])){
      transcrip_edge_num <- transcrip_edge_num+1
      transcrip_adjacency$Node1[transcrip_edge_num] <- i-1
      transcrip_adjacency$Node2[transcrip_edge_num] <- j-1
    }
  }
}

final_transcrip_adjacency <- data.frame(Node1=transcrip_adjacency$Node1, Node2=transcrip_adjacency$Node2)
saveRDS(final_transcrip_adjacency, file=paste0(output_path, "transcrip_adjacency.rds"))
write.csv(final_transcrip_adjacency, paste0(output_path, "transcrip_adjacency.csv"), row.names=FALSE)
