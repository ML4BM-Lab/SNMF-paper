
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
output_path <- args[1]
raw_proportions_path <- args[2]
proportions_path <- args[3]

raw_proportions <- t(read.csv(raw_proportions_path, row.names=1))
proportions <- read.csv(proportions_path, row.names=1)

library(RcppHungarian)

k <- dim(raw_proportions)[1]

cost <- as.matrix(dist(rbind(raw_proportions, t(proportions))))[(k+1):(2*k),1:k]
hungarian_result <- HungarianSolver(cost)
I <- hungarian_result$pairs[,2]
raw_proportions <- t(raw_proportions[I,])
rownames(raw_proportions) <- rownames(proportions)
colnames(raw_proportions) <- colnames(proportions)

write.csv(raw_proportions, paste0(output_path, "hungarian_proportions.csv"))