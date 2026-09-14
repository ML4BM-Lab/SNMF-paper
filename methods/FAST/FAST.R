
library(FAST)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
output_path <- args[1]
seed <- as.integer(args[2])

set.seed(seed)

load(paste0(output_path, "tmp/X.RData"))
load(paste0(output_path, "tmp/adj.RData"))
load(paste0(output_path, "tmp/config.RData"))

res <- dmain(X, adj, config)
H <- res$H
rownames(H) <- colnames(X)
save(H, file=paste0(output_path, "tmp/raw_proportions.RData"))
