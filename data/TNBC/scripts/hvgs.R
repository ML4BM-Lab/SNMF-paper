script_arg <- grep("--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("--file=", "", script_arg[1])))
} else {
  getwd()
}
data_dir <- normalizePath(file.path(script_dir, ".."))
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."))
zenodo_dir <- file.path(repo_root, "data", "zenodo", "TNBC")
dir.create(zenodo_dir, recursive = TRUE, showWarnings = FALSE)

# counts: genes x spots matrix
counts <- read.csv(file.path(zenodo_dir, "TNBC_counts.csv"), row.names=1, check.names=FALSE)

N <- 5000

gene_var <- apply(counts, 1, var)
top <- names(sort(gene_var, decreasing = TRUE))[1:N]

markers <- read.csv(file.path(zenodo_dir, "TNBC_marker_genes.csv"))
marker_genes <- unique(markers$gene)

genes_keep <- union(top, marker_genes)

counts_subset <- counts[genes_keep, ]

cat(
  "HVGs:", length(top),
  "\nMarkers found:", length(marker_genes),
  "\nGenes retained:", length(genes_keep), "\n"
)

write.csv(counts_subset, file.path(zenodo_dir, paste0("TNBC_counts_hvgs", N, ".csv")))
