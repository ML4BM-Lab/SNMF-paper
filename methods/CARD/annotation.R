# For the interpretation, one of the exploratory analysis we have done to link the cell types of the marker genes and the CT labels is as the following 
# (You can also see supplementary figure 83 in the paper):

#  For each "CT" label, i.e. CT14 here in the supplementary Figure 83, we first divide the spots into CT14 enriched and CT14 non-enriched locations, 
# and the CT14 enriched spatial locations were defined as the spatial locations that contains at least 50% of CT14 cell type (You can also define it as 
# the median or > 50%). And then for each set of cell type specific marker genes, we calculate the mean gene expression of each set of marker genes in 
# the CT14 enriched locations, the cell type corresponds to the highest set of mean gene expression can be assumed as the cell type annotation for the 
# CT14 label. And then we iterate the above procedure for each "CT" label.

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
data_path <- args[1]
marker_ref_path <- args[2]
output_path <- args[3]

data <- read.csv(data_path, row.names=1, check.names = FALSE)
data <- sweep(data, 2, colSums(data), FUN = "/") * 1e4 # Normalization

marker_ref_df <- read.csv(marker_ref_path)
marker_list <- split(marker_ref_df$gene, marker_ref_df$cluster)

load(paste0(output_path, "tmp/raw_proportions.RData"))

# Option 1: CT_* choose cell types
# mapping <- list()
# for (i in colnames(proportions)) {
#     threshold <- 0.5
#     enriched_spots <- rownames(proportions[proportions[, i] > threshold,])

#     # gradually decrease threshold until we have enriched spots
#     while (length(enriched_spots) == 0 && threshold > 0) {
#         threshold <- threshold - 0.05
#         enriched_spots <- rownames(proportions[proportions[, i] > threshold,])
#     }

#     max_exp <- -1
#     max_ct <- NA
#     for (ct in names(marker_list)) {
#         if (ct %in% as.vector(unlist(mapping))) {
#             next
#         }
#         gene_set <- marker_list[[ct]]
#         expression <- data[gene_set[gene_set %in% rownames(data)], enriched_spots]
#         avg_expression <- sum(expression) / sum(gene_set %in% rownames(data))
#         if (avg_expression > max_exp) {
#             max_exp <- avg_expression
#             max_ct <- ct
#         }
#     }

#     mapping[i] <- max_ct
# }

# Option 2: cell_types choose CT_*
mapping <- list()
for (ct in names(marker_list)) {
    max_exp <- -1
    max_comp <- NA

    for (i in colnames(proportions)) {
        if (i %in% names(mapping)) {
            next
        }
 
        threshold <- 0.5
        enriched_spots <- rownames(proportions[proportions[, i] > threshold,])

        # gradually decrease threshold until we have enriched spots
        while (length(enriched_spots) == 0 && threshold > 0) {
            threshold <- threshold - 0.05
            enriched_spots <- rownames(proportions[proportions[, i] > threshold,])
        }

        gene_set <- marker_list[[ct]]
        expression <- data[gene_set[gene_set %in% rownames(data)], enriched_spots, drop = FALSE]
        avg_expression <- sum(expression) / sum(gene_set %in% rownames(data))

        if (avg_expression > max_exp) {
            max_exp <- avg_expression
            max_comp <- i
        }
    }
    
    mapping[max_comp] <- ct
}

colnames(proportions) <- mapping[colnames(proportions)]

write.csv(proportions, paste0(output_path, "CARD_proportions.csv"))
