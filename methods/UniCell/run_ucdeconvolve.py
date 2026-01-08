
import sys
import os

import ucdeconvolve as ucd
import scanpy as sc

output_path = sys.argv[1]
# seed = sys.argv[2]

adata = sc.read_h5ad(os.path.join(output_path, "tmp", "adata.h5ad"))

token = "uc_MbArWLfkpieqtYv1fIQbNMlKuKCW4IvfFDzKt35t03mNRpaS"
ucd.api.authenticate(token)
ucd.tl.base(adata)

predictions = ucd.utils.read_results(adata, category = 'raw')
predictions.to_csv(os.path.join(output_path, "raw_predictions.csv"))