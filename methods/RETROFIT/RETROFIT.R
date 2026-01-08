
library(retrofit)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
output_path <- args[1]
seed <- as.integer(args[2])

set.seed(seed)

load(paste0(output_path, "tmp/data.RData"))
res <- retrofit::decompose(x, L=32, iterations=100, verbose=TRUE)
save(res, file=paste0(output_path, "tmp/result.RData"))
