
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
data_path <- args[1]
marker_ref_path <- args[2]
output_path <- args[3]

# Load data
x <- read.csv(data_path, row.names=1, check.names=FALSE)
x <- t(x)
colnames(x) <- paste0(toupper(substr(colnames(x), 1, 1)), tolower(substr(colnames(x), 2, nchar(colnames(x)))))
save(x, file=paste0(output_path, "tmp/data.RData"))

marker_ref_df <- read.csv(marker_ref_path)
marker_ref <- split(marker_ref_df$gene, marker_ref_df$cluster)
marker_ref <- lapply(marker_ref, function(x) {
  sapply(x, function(g) {
    g_lower <- tolower(g)               # all lowercase
    paste0(toupper(substring(g_lower,1,1)), substring(g_lower,2))  # capitalize first
  })
})

overlap <- intersect(colnames(x), unique(unlist(marker_ref)))
marker_counts <- colSums(x[,overlap])
low_expr_markers <- overlap[marker_counts == 0]
keeping_genes <- setdiff(overlap, low_expr_markers)
marker_ref <- lapply(marker_ref, function(gs) {
  gs <- intersect(gs, keeping_genes)   # keep only valid, expressed genes
  gs                                   # return cleaned list
})

save(marker_ref, file=paste0(output_path, "tmp/markers.RData"))