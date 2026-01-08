
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
marker_ref_path <- args[1]

marker_ref_df <- read.csv(marker_ref_path, row.names=1)
marker_ref <- lapply(1:nrow(marker_ref_df), function(i) as.vector(unlist(marker_ref_df[i, ])))
names(marker_ref) <- rownames(marker_ref_df)

load(paste0(output_path, "tmp/proportions.RData"))

celltype_annotations <- annotateCellTypesGSEA(beta = results$beta, gset = marker_ref, qval = 0.05)
proportions <- results$theta
colnames(proportions) <- celltype_annotations[colnames(proportions)]