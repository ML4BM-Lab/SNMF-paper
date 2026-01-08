
import sys
import os

import pandas as pd
import scanpy as sc

data_path = sys.argv[1]
output_path = sys.argv[2]

## Load data
data = pd.read_csv(data_path, index_col=0)
adata = sc.AnnData(
    X=data.values.T,
    obs=pd.DataFrame(index=data.columns),  # rows are observations
    var=pd.DataFrame(index=data.index)     # rows are variables
)
adata.obs[["array_row", "array_col"]] = adata.obs_names.str.split("x").tolist()

adata.write_h5ad(os.path.join(output_path, "tmp", "adata.h5ad"))