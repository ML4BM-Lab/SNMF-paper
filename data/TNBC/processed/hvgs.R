# counts: genes x spots matrix

counts <- read.csv("/scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/counts.csv", row.names=1, check.names=FALSE)

N <- 5000

gene_var <- apply(counts, 1, var)
top <- names(sort(gene_var, decreasing = TRUE))[1:N]

markers <- read.csv("/scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/marker_genes.csv")
marker_genes <- unique(markers$gene)

genes_keep <- union(top, marker_genes)

counts_subset <- counts[genes_keep, ]

cat(
  "HVGs:", length(top),
  "\nMarkers found:", length(marker_genes),
  "\nGenes retained:", length(genes_keep), "\n"
)

write.csv(counts_hvg, paste0("/scratch/lalonsoeste/PhD/NMF_deconvolution/data/TNBC/processed/counts_hvgs", N, ".csv"))