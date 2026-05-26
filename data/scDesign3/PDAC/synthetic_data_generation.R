
library(scDesign3)
library(SingleCellExperiment)
library(scran)
library(ggplot2)
library(dplyr)
library(viridis)
library(IOBR)
library(scatterpie)
theme_set(theme_bw())

load("./sc_count.RData")
load("./sc_meta.RData")

sc_sce <- SingleCellExperiment(assays = list(counts = sc_count), colData = sc_meta)

dec <- modelGeneVar(sc_sce, assay.type = "counts")
top_hvg <- getTopHVGs(dec, n = 1000)
sc_sce <- sc_sce[top_hvg, ]

load("./spatial_count.RData")
load("./spatial_location.RData")

spatial_sce <- SingleCellExperiment(assays = list(counts = spatial_count), colData = spatial_location)
spatial_sce <- spatial_sce[top_hvg[top_hvg %in% rownames(spatial_sce)],]

cell_type <- unique(colData(sc_sce)$cellType)

# Simulation
set.seed(123)

# Cell-type reference estimation
sc_data <- construct_data(
  sce = sc_sce,
  assay_use = "counts",
  celltype = "cellType",
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = NULL,
  corr_by = "1"
)

sc_marginal <- fit_marginal(
  data = sc_data,
  predictor = "gene",
  mu_formula = "cellType",
  sigma_formula = "cellType",
  family_use = "nb",
  n_cores = 1,
  usebam = FALSE,
  parallelization = "pbmcmapply"
  
)

sc_copula <- fit_copula(
  sce = sc_sce,
  assay_use = "counts",
  marginal_list = sc_marginal,
  family_use = "nb",
  copula = "gaussian",
  n_cores = 1,
  input_data = sc_data$dat
)

sc_para <- extract_para(
  sce = sc_sce,
  marginal_list = sc_marginal,
  n_cores = 1,
  family_use = "nb",
  new_covariate = sc_data$newCovariate,
  data = sc_data$dat
)

sc_newcount <- simu_new(
  sce = sc_sce,
  mean_mat = sc_para$mean_mat,
  sigma_mat = sc_para$sigma_mat,
  zero_mat = sc_para$zero_mat,
  quantile_mat = NULL,
  copula_list = sc_copula$copula_list,
  n_cores = 1,
  family_use = "nb",
  input_data = sc_data$dat,
  new_covariate = sc_data$newCovariate,
  filtered_gene = sc_data$filtered_gene
)



set.seed(123)
spatial_data <- construct_data(
  sce = spatial_sce,
  assay_use = "counts",
  celltype = NULL,
  pseudotime = NULL,
  spatial = c("x", "y"),
  other_covariates = NULL,
  corr_by = "1"
)

spatial_marginal <- fit_marginal(
  data = spatial_data,
  predictor = "gene",
  mu_formula = "s(x, y, bs = 'gp', k = 50, m = c(1, 2, 1))",
  sigma_formula = "1",
  family_use = "nb",
  n_cores = 2,
  usebam = FALSE, 
  parallelization = "pbmcmapply"
)

spatial_copula <- fit_copula(
  sce = spatial_sce,
  assay_use = "counts",
  marginal_list = spatial_marginal,
  family_use = "nb",
  copula = "gaussian",
  n_cores = 1,
  input_data = spatial_data$dat
)

spatial_para <- extract_para(
  sce = spatial_sce,
  marginal_list = spatial_marginal,
  n_cores = 2,
  family_use = "nb",
  new_covariate = spatial_data$newCovariate,
  data = spatial_data$dat
)


#Now we get the fitted models for scRNA-seq and spatial data. We need to extract their mean parameters (i.e., expected expression values).
sc_sig_matrix <- sapply(cell_type, function(x) {
  rowMeans(t(sc_para$mean_mat)[, colData(sc_sce)$cellType %in% x])
})

spatial_matrix <- (t(spatial_para$mean_mat))

# We use CIBERSORT to decompose each spot’s expected expression into cell-type proportions. This step is to set the true cell-type proportions.
sig_matrix <- as.data.frame(sc_sig_matrix)
mixture_file <- as.data.frame(spatial_matrix)

proportion_mat <- IOBR::CIBERSORT(sig_matrix, mixture_file, QN = FALSE, absolute = FALSE, perm = 10)
proportion_mat <- proportion_mat[,1:(ncol(proportion_mat)-3)]

write.csv(proportion_mat, "./counts/proportions.csv", row.names=T)

# Pie chart visualization
colors_cell_type <- c("#E69F00", "#56B4E9", "#009E73", 
                      "#0072B2")
d_pie <- as_tibble(colData(spatial_sce), rownames = "cell") %>% bind_cols(as_tibble(proportion_mat)) %>% dplyr::mutate(region = seq_len(dim(spatial_sce)[2])) %>% dplyr::mutate(X= x, Y = y)

p_pie_plot <- ggplot() + geom_scatterpie(aes(x=X, y=Y, group=region), data=d_pie ,
                                         cols = cell_type, color=NA) + coord_fixed(ratio = 1) + 
  scale_fill_manual(values = colors_cell_type) + coord_equal()+ theme_bw() + theme(legend.position = "left")  + theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())+ guides(fill=guide_legend(title="Cell type"))
p_pie_plot


# Then we can simulate new spatial data where each spot is the sum of 50 cells/5 (therefore on average 10 cells per spot). Increasing the number of cells will make the spatial data smoother (closer to the expected spatial expression).
set.seed(123)
scsim_sce <- sc_sce
counts(scsim_sce) <- sc_newcount

spatial_new_mixture <- (apply(proportion_mat, 1, function(x) {
  n = 50
  rowSums(sapply(cell_type, function(y) {
    index <- sample(which(colData(scsim_sce)$cellType==y), size = n, replace = length(which(colData(scsim_sce)$cellType==y)) < 50)
    rowSums(sc_newcount[, index]) * x[y]
  }))
}))

spatial_new_mixture <- spatial_new_mixture/5

### Ceiling to integer
spatial_new_mixture <- ceiling(spatial_new_mixture)

write.csv(spatial_new_mixture, "./counts/mixture_file.csv", row.names=T)

spatialmix_sce <- spatial_sce
counts(spatialmix_sce) <- as.matrix(spatial_new_mixture)


#Finally, we can check the simulated results. We use four cell-type marker genes as the example.
sc_sig_matrix <- sapply(cell_type, function(x) {
  rowMeans(t(sc_para$mean_mat)[, colData(sc_sce)$cellType %in% x])
})
spatial_sc_mixture <- tcrossprod(as.matrix(sc_sig_matrix), as.matrix(proportion_mat))

rownames(spatial_sc_mixture) <- rownames(spatial_new_mixture)

location <- colData(spatial_sce)
spatial_real_tbl <- as_tibble(t(log1p(counts(spatial_sce)))) %>% dplyr::mutate(X = location$x,
                                                                               Y = location$y) %>%
  tidyr::pivot_longer(-c("X", "Y"), names_to = "Gene", values_to = "Expression") %>% dplyr::mutate(Method = "Real data")

spatial_real_tbl <- transform(spatial_real_tbl, Expression=ave(Expression, Gene, FUN=scales::rescale))


spatial_mixture_tbl <- as_tibble(t(log1p(spatial_new_mixture))) %>% dplyr::mutate(X = location$spatial1,
                                                                                  Y = location$spatial2) %>%
  tidyr::pivot_longer(-c("X", "Y"), names_to = "Gene", values_to = "Expression") %>% dplyr::mutate(Method = "scDesign3")

spatial_mixture_tbl <- transform(spatial_mixture_tbl, Expression=ave(Expression, Gene, FUN=scales::rescale))

spatial_tbl <- bind_rows(list(spatial_real_tbl, spatial_mixture_tbl))

sc_marker <- c("Penk", "Apold1", "Cdhr1", "S100a5")

p_prop <- spatial_tbl %>% dplyr::filter(Gene %in% sc_marker) %>% dplyr::mutate(Gene = factor(Gene, levels = sc_marker)) %>% ggplot(aes(x = X, y = Y, color = Expression))  + ggrastr::rasterize(geom_point(size = 1), dpi = 300) + scale_colour_gradientn(colors = viridis_pal(option = "B", direction = -1)(10), limits=c(0, 1)) + coord_fixed(ratio = 1) + facet_grid(Method ~ Gene ) + theme_bw() + theme(legend.position = "right")  + theme(
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks.y=element_blank())
p_prop
