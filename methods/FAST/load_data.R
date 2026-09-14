
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
data_path <- args[1]
output_path <- args[2]
k <- as.integer(args[3])

X <- read.csv(data_path, row.names=1, check.names=FALSE)

X <- log(X+1)
X <- X/max(X)

# adj
spots <- colnames(X)

# Parse coordinates
coords <- do.call(rbind, strsplit(spots, "x"))
coords <- apply(coords, 2, as.integer)

# Create spots x spots matrix
adj <- matrix(
  0,
  nrow = length(spots),
  ncol = length(spots),
  dimnames = list(spots, spots)
)

# Fill adjacency matrix
for (i in seq_along(spots)) {
  for (j in seq_along(spots)) {
    
    # Manhattan distance
    d <- sum(abs(coords[i, ] - coords[j, ]))
    
    if (d == 1) {
      adj[i, j] <- 1
    }
  }
}

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