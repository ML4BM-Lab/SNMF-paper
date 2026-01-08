
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
output_path <- args[1]
k <- as.integer(args[2])
num_initializations <- as.integer(args[3])
probs <- as.numeric(args[4])
seed <- as.integer(args[5])

library(GPUmatrix)
library(torch)

load(paste0(output_path, "tmp/data.RData"))
load(paste0(output_path, "tmp/S.RData"))

source("/scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SNMF/NMFKLMixing.R")

counts <- gpu.matrix(counts, dtype = "float32")
S <- gpu.matrix(S, dtype = "float32")

set.seed(seed)

output <- NMFKLMixing(counts, S = S, k = k,
                      niter = 2000, tol = 1e-4, num_initializations=num_initializations)

W <- as.matrix(output$W)
H <- as.matrix(output$H)

D <- diag(matrixStats::colQuantiles(W, probs = probs, na.rm = T))
D_1 <- diag(1/matrixStats::colQuantiles(W, probs = probs, na.rm = T))

# Normalize the W and H matrices
W <- W %*% D_1; 
H <- D %*% H
HC <- as.matrix(H %*% S)

# To get an H matrix of proportions.
HC <- t(t(HC)/colSums(HC))

save(W, file=paste0(output_path, "tmp/W.RData"))
save(HC, file=paste0(output_path, "tmp/raw_proportions.RData"))
