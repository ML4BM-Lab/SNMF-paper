
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
data_path <- args[1]
output_path <- args[2]
gamma <- as.numeric(args[3])

counts <- read.csv(data_path, row.names=1, check.names=FALSE)

filter <- rowSums(counts) > 10
counts <- counts[filter, ]
counts <- as.matrix(counts)

positions <- matrix(as.numeric(unlist(strsplit(colnames(counts), "x"))), 
                ncol=2, byrow = TRUE)
x <- positions[,1]
y <- positions[,2]

S <- exp(-gamma * as.matrix(dist(cbind(x,y)))^2)
S[S < 1e-3] <- 0 

S <- Diagonal(x = 1/rowSums(S)) %*% S

save(counts, file=paste0(output_path, "tmp/data.RData"))
save(S, file=paste0(output_path, "tmp/S.RData"))