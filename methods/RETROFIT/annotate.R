
library(retrofit)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
output_path <- args[1]
marker_ref_path <- args[2]

load(paste0(output_path, "tmp/result.RData"))

W <- res$w
H <- res$h

## match the latent components to known cell types using a list of cell-type-specific marker genes
marker_ref_df <- read.csv(marker_ref_path)
marker_ref <- split(marker_ref_df$gene, marker_ref_df$cluster)

res <- retrofit::annotateWithMarkers(marker_ref, K=length(marker_ref), decomp_w=W, decomp_h=H)
H_annotated <- res$h
proportions <- sweep(H_annotated, 2, colSums(H_annotated), FUN = "/")

write.csv(t(proportions), paste0(output_path, "RETROFIT_proportions.csv"))