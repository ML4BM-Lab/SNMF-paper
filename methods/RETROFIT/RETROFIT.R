
library(retrofit)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
output_path <- args[1]
k <- as.integer(args[2])
seed <- as.integer(args[3])

set.seed(seed)

load(paste0(output_path, "tmp/data.RData"))
res <- retrofit::decompose(x, L=2*k, iterations=100, verbose=TRUE)
save(res, file=paste0(output_path, "tmp/result.RData"))
