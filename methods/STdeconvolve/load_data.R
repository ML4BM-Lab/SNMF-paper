
library(STdeconvolve)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
data_path <- args[1]
output_path <- args[2]

# Load data
x <- as.matrix(read.csv(data_path, row.names=1, check.names=FALSE))

counts <- cleanCounts(x, min.lib.size = 100)
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)

save(corpus, file=paste0(output_path, "tmp/corpus.RData"))