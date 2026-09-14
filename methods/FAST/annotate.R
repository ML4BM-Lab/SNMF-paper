
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
output_path <- args[1]
proportions_path <- args[2]

proportions <- read.csv(proportions_path, row.names=1)
load(paste0(output_path, "tmp/raw_proportions.RData"))

library(RcppHungarian)

k <- dim(H)[2]

cost <- as.matrix(dist(rbind(t(H), t(proportions))))[(k+1):(2*k),1:k]
hungarian_result <- HungarianSolver(cost)
I <- hungarian_result$pairs[,2]
H <- H[,I]
colnames(H) <- colnames(proportions)

write.csv(H, paste0(output_path, "FAST_proportions.csv"))