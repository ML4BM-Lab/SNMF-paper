
library(CARD)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
output_path <- args[1]

load(paste0(output_path, "tmp/CARD_obj.RData"))
CARDfree_obj = CARD_refFree(CARDfree_obj)
write.csv(CARDfree_obj@Proportion_CARD, paste0(output_path, "CARD_proportions.csv"))













