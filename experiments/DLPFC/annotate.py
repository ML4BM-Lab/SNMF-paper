
import sys

import scanpy as sc
import pandas as pd

data_path = sys.argv[1]
proportions_path = sys.argv[2]

adata = sc.read_h5ad(data_path)

proportions = pd.read_csv(proportions_path)
