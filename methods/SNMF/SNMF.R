
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
output_path <- args[1]
loss_func <- args[2]
dispersion_mode <- args[3]
k <- as.integer(args[4])
seed <- as.integer(args[5])

library(GPUmatrix)
library(torch)

load(paste0(output_path, "tmp/data.RData"))
load(paste0(output_path, "tmp/S.RData"))

if (loss_func == "KL_poisson") {
    source("/scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SNMF/update_rules/KL_poisson.R")
} else if (loss_func == "squared_error") {
    source("/scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SNMF/update_rules/squared_error.R")
} else if (loss_func == "KL_NB") {
    source("/scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SNMF/update_rules/KL_NB.R")
} else {
    stop("Loss function selected not supported.")
}

gpu_counts <- gpu.matrix(counts, dtype = "float32")
S <- gpu.matrix(S, dtype = "float32")

set.seed(seed)

output <- snmf(gpu_counts, S = S, k = k,
                      niter = 2000, tol = 1e-4, num_initializations=10, dispersion_mode=dispersion_mode)

W <- as.matrix(output$W)
H <- as.matrix(output$H)

D <- diag(matrixStats::colQuantiles(W, probs = 0.75, na.rm = T))
D_1 <- diag(1/matrixStats::colQuantiles(W, probs = 0.75, na.rm = T))

# Normalize the W and H matrices
W <- W %*% D_1; 
H <- D %*% H
HC <- as.matrix(H %*% S)

# To get an H matrix of proportions.
HC <- t(t(HC)/colSums(HC))

save(W, file=paste0(output_path, "tmp/W.RData"))
save(HC, file=paste0(output_path, "tmp/raw_proportions.RData"))

# Name spots and genes
colnames(HC) <- colnames(counts)
HC <- t(HC)

rownames(W) <- rownames(counts)

write.csv(W, paste0(output_path, "signatures.csv"))
write.csv(HC, paste0(output_path, "SNMF_proportions.csv"))

write(output$niter, file=paste0(output_path, "niter.txt"))

if (loss_func == "KL_NB") {
    phi <- as.matrix(output$phi)

    rownames(phi) <- rownames(counts)
    colnames(phi) <- colnames(counts)

    write.csv(phi, paste0(output_path, "phi.csv"))

    alpha <- matrix(output$alpha, nrow = 1)
    beta <- matrix(output$beta, nrow = 1)

    colnames(alpha) <- rownames(counts)
    colnames(beta) <- colnames(counts)

    write.csv(alpha, paste0(output_path, "alpha.csv"))
    write.csv(beta, paste0(output_path, "beta.csv"))
}
