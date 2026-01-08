
library(SMART)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
output_path <- args[1]
seed <- as.integer(args[2])

load(paste0(output_path, "tmp/data.RData"))
load(paste0(output_path, "tmp/markers.RData"))

base_res <- SMART_base(stData=x,
                       markerGs=marker_ref, 
                       noMarkerCts=0, outDir=output_path, seed=seed)

write.csv(base_res$ct_proportions, paste0(output_path, "SMART_proportions.csv"))