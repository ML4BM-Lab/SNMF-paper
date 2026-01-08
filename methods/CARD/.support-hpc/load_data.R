
library(CARD)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
data_path <- args[1]
marker_ref_path <- args[2]
output_path <- args[3]

data <- read.csv(data_path, row.names=1, check.names = FALSE)

spatial_count <- as.matrix(data)
spatial_location <- do.call(rbind, strsplit(colnames(data), "x"))
spatial_location <- data.frame(x = as.numeric(spatial_location[,1]), y = as.numeric(spatial_location[,2]))
rownames(spatial_location) <- colnames(spatial_count)

marker_list_df <- read.csv(marker_ref_path, row.names=1)
marker_list <- lapply(1:nrow(marker_list_df), function(i) as.vector(unlist(marker_list_df[i, ])))
names(marker_list) <- rownames(marker_list_df)

CARDfree_obj = createCARDfreeObject(
  markerList = marker_list,
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  minCountGene = 1,
  minCountSpot = 1)

save(CARDfree_obj, file=paste0(output_path, "tmp/CARD_obj.RData"))