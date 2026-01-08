
library(CARD)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
output_path <- args[1]
seed <- as.integer(args[2])

set.seed(seed)

load(paste0(output_path, "tmp/CARD_obj.RData"))
CARDfree_obj = CARD_refFree(CARDfree_obj)
proportions <- CARDfree_obj@Proportion_CARD
save(proportions, file=paste0(output_path, "tmp/raw_proportions.RData"))
