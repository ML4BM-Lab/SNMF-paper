
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) stop("Please provide the data directory path")
data_path <- args[1]
output_path <- args[2]

# Load data
x <- read.csv(data_path, row.names=1, check.names=FALSE)
save(x, file=paste0(output_path, "tmp/data.RData"))

print("Data loaded and saved")