# --- Libraries -------------------------------------------------------------
library(tidyverse)
library(pheatmap)
library(UpSetR)
library(RColorBrewer)
library(ggrepel)

# --- Paths -----------------------------------------------------------------
marker_genes_path <- "/scratch/lalonsoeste/PhD/NMF_deconvolution/data/scDesign3/PDAC/marker_genes/sig_marker_genes.csv"
signatures_path   <- "/scratch/lalonsoeste/PhD/NMF_deconvolution/results/benchmarking/PDAC/SNMF/signatures.csv"

# Create output directory
out_dir <- "/scratch/lalonsoeste/PhD/NMF_deconvolution/results/benchmarking/PDAC/SNMF/signatures_analysis/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load data -------------------------------------------------------------
marker_ref_df <- read.csv(marker_genes_path)
marker_ref <- split(marker_ref_df$gene, marker_ref_df$cluster)

n_markers <- min(lengths(marker_ref))   # equalize to smallest list size
marker_ref <- lapply(marker_ref, function(genes) head(genes, n_markers))

signatures <- read.csv(signatures_path, row.names = 1, check.names = FALSE)
signatures_mat <- as.matrix(signatures)

# --- 2. Dotplot of top genes per signature --------------------------------
top_genes_per_sig <- apply(signatures_mat, 2, function(x) {
  head(names(sort(x, decreasing = TRUE)), 3)
}) %>%
  unlist() %>%
  unique()

plot_df <- signatures_mat[top_genes_per_sig, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "Signature", values_to = "Weight")

p_dot <- ggplot(plot_df, aes(Signature, gene, size = Weight, color = Weight)) +
  geom_point() +
  scale_color_viridis_c() +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  # <- rotate x labels
  ) +
  ggtitle("Top weighted genes per signature")

ggsave(file.path(out_dir, "top_genes_dotplot.pdf"), p_dot, width = 8, height = 10)

# --- 3. UpSet plot ---------------------------------------------------------
top_n <- 30

# Get top N genes per signature based on their weights
top_genes_signatures <- lapply(seq_len(ncol(signatures_mat)), function(i) {
  x <- signatures_mat[, i]
  # ensure we keep gene names (rownames of signatures_mat)
  genes_sorted <- names(sort(x, decreasing = TRUE))
  genes_sorted[seq_len(min(top_n, length(genes_sorted)))]
})
names(top_genes_signatures) <- colnames(signatures_mat)

# Build binary matrix: genes × signatures (TRUE if gene in top N for that signature)
all_top_genes <- unique(unlist(top_genes_signatures))
sig_binary <- sapply(top_genes_signatures, function(genes) as.integer(all_top_genes %in% genes))
sig_binary <- as.data.frame(sig_binary)
rownames(sig_binary) <- all_top_genes

pdf(file.path(out_dir, "marker_overlap_upset.pdf"), width = 9, height = 6)
upset(sig_binary,
      sets = colnames(sig_binary),
      order.by = "freq",
      mainbar.y.label = "Gene Overlaps Among Signatures (Top 30 Genes)",
      sets.x.label = "Genes per Signature (Top 30)")
dev.off()

# --- 4. Correlation heatmap between marker sets and signatures ------------
marker_summary <- lapply(marker_ref, function(genes) intersect(genes, rownames(signatures_mat)))
marker_sig_scores <- sapply(marker_summary, function(genes) {
  if (length(genes) == 0) return(rep(NA, ncol(signatures_mat)))
  colMeans(signatures_mat[genes, , drop = FALSE])
}) # Marker sets x Signatures (5x5)

corr_mat <- cor(marker_sig_scores, use = "pairwise.complete.obs", method = "pearson") # Pairwise correlation of marker sets

pdf(file.path(out_dir, "marker_signature_correlation.pdf"), width = 10, height = 8)
pheatmap(corr_mat, color = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
         main = "Correlation between Marker Sets and Signatures")
dev.off()

# --- 5. PCA visualization --------------------------------------------------
pca <- prcomp(t(signatures_mat), scale. = TRUE)
pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("Signature")

marker_means <- sapply(marker_ref, function(genes) {
  present <- intersect(genes, rownames(signatures_mat))
  colMeans(signatures_mat[present, , drop = FALSE])
})

marker_means_df <- as.data.frame(marker_means)
rownames(marker_means_df) <- colnames(signatures_mat)

pca_df$Dominant_marker <- colnames(marker_means_df)[apply(marker_means_df, 1, which.max)]

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Dominant_marker, label = Signature)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3.5, max.overlaps = 15) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  ggtitle("PCA of signatures (colored by dominant marker set)")

ggsave(file.path(out_dir, "signature_pca.pdf"), p_pca, width = 8, height = 6)

# --- Summary message ------------------------------------------------------
cat("\n✅ All figures saved in:\n", out_dir, "\n")
