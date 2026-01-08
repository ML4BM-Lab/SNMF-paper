
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
output_path <- args[1]
proportions_path <- args[2]

proportions <- read.csv(proportions_path, row.names=1)

load(paste0(output_path, "tmp/data.RData"))
load(paste0(output_path, "tmp/W.RData"))
load(paste0(output_path, "tmp/raw_proportions.RData"))

library(RcppHungarian)

k <- dim(HC)[1]

cost <- as.matrix(dist(rbind(HC, t(proportions))))[(k+1):(2*k),1:k]
hungarian_result <- HungarianSolver(cost)
I <- hungarian_result$pairs[,2]
HC <- t(HC[I,])
rownames(HC) <- rownames(proportions)
colnames(HC) <- colnames(proportions)

W <- W[,I]
rownames(W) <- rownames(counts)
colnames(W) <- colnames(proportions)

write.csv(HC, paste0(output_path, "SNMF_proportions.csv"))
write.csv(W, paste0(output_path, "signatures.csv"))