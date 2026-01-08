# Recovered: not tested yet. Something might fail

library(NMF)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript nmf.R <data.csv> <output_dir/> <k> <nrun>")
}

data_path  <- args[1]
output_path <- args[2]
k     <- as.integer(args[3])
nrun  <- as.integer(args[4])

# Ensure output directory ends with /
if (!grepl("/$", output_path)) {
  output_path <- paste0(output_path, "/")
}

# Load data
data <- read.csv(data_path, row.names = 1, check.names = FALSE)

# Run NMF with KL divergence
result <- nmf(
  data,
  rank = k,
  method = "KL",
  nrun = nrun,
  seed = 42
)

# Extract matrices
W <- basis(result)   # features x k
H <- coef(result)    # k x samples

# Save results
write.csv(W, file = paste0(output_path, "W.csv"))
write.csv(H, file = paste0(output_path, "H.csv"))
