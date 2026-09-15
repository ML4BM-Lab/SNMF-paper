
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
data_path <- args[1]
output_path <- args[2]
k <- as.integer(args[3])

X <- read.csv(data_path, row.names=1, check.names=FALSE)
X <- X[rowSums(X) > 0, , drop = FALSE]

X <- log(X+1)
X <- X/max(X)

# adj
spots <- colnames(X)

# Parse coordinates
coords <- do.call(rbind, strsplit(spots, "x"))
coords <- apply(coords, 2, as.integer)

# Create spots x spots matrix
distance <- matrix(
  0,
  nrow = length(spots),
  ncol = length(spots),
  dimnames = list(spots, spots)
)

for (i in seq_along(spots)) {
  for (j in seq_along(spots)) {
    distance[i, j] <- (coords[i,1] - coords[j,1])^2 + (coords[i,2] - coords[j,2])^2
  }
}

adj <- matrix(
  0,
  nrow = length(spots),
  ncol = length(spots),
  dimnames = list(spots, spots)
)

for (i in seq_along(spots)) {
  d <- distance[i, ]
  nn <- order(d)[seq_len(5)]

  adj[i, nn] <- 1
}

# Make adj symmetric
adj <- pmax(adj, t(adj))

config <- list(
    r = k,
    lambda_1 = 0.05,
    lambda_2 = 0.1,
    n_iter = 1000,
    converge = 2e-06
)

save(X, file=paste0(output_path, "tmp/X.RData"))
save(adj, file=paste0(output_path, "tmp/adj.RData"))
save(config, file=paste0(output_path, "tmp/config.RData"))