
library(STdeconvolve)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
output_path <- args[1]
marker_ref_path <- args[2]
minK <- as.integer(args[3])
hungarian <- args[4] == "true"
seed <- as.integer(args[5])

# Load data
load(paste0(output_path, "tmp/corpus.RData"))

marker_ref_df <- read.csv(marker_ref_path)
marker_ref <- split(marker_ref_df$gene, marker_ref_df$cluster)

if (hungarian) {
    ldas <- fitLDA(t(as.matrix(corpus)), Ks = minK, plot=FALSE, verbose=FALSE, seed=seed)
    optLDA <- optimalModel(models = ldas, opt = minK)
    results <- getBetaTheta(optLDA, perc.filt = 0.0, betaScale = 1000)
} else {
    K <- minK
    n_cell_types <- 0
    while (n_cell_types < minK - 1) {
        if (K > 100) {
            ldas <- fitLDA(t(as.matrix(corpus)), Ks = K, plot=FALSE, verbose=FALSE, seed=seed, ncores=parallel::detectCores())
            optLDA <- optimalModel(models = ldas, opt = K)
            results <- getBetaTheta(optLDA, perc.filt = 0.0, betaScale = 1000)
            break
        } else {
            print(paste0("Running STdeconvolve for K=", K))
        }

        ## choose optimal number of cell-types
        ldas <- fitLDA(t(as.matrix(corpus)), Ks = K, plot=FALSE, verbose=FALSE)
        ## get best model results
        optLDA <- optimalModel(models = ldas, opt = K)

        ## extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
        results <- getBetaTheta(optLDA, perc.filt = 0.0, betaScale = 1000)

        celltype_annotations <- annotateCellTypesGSEA(beta = results$beta, gset = marker_ref, qval = 0.05)
        results <- celltype_annotations$results

        mapping <- list()
        for (k in seq(1,K)) {
            annot <- results[[as.character(k)]]
            if (dim(annot)[1]) {
                mapping[[as.character(k)]] <- as.character(annot[order(annot[["padj"]])][1,"pathway"])
            } else {
                mapping[[as.character(k)]] <- NA
            } 
        }
        n_cell_types <- length(unique(mapping[!is.na(mapping)]))

        print(paste0("Detected ", n_cell_types, " cell types"))

        K <- K + 1
    }
}

write.csv(results$theta, paste0(output_path, "STdeconvolve_proportions.csv"))
