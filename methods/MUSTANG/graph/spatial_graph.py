#!/usr/bin/env python3
import sys
import pandas as pd
from bayestme import utils, data

data_path = sys.argv[1]
output_path = sys.argv[2]
technology = sys.argv[3]

if ".csv" in data_path:
    pos = pd.DataFrame(pd.read_csv(data_path, index_col=0).columns.str.split("x").tolist(), columns=['x','y'])
neighbors = utils.get_edges(pos.values, layout=data.Layout.HEX if technology == "Visium" else data.Layout.SQUARE)

neighbors_file = f"{output_path}/spatial_neighbors.csv"
pd.DataFrame(neighbors).to_csv(neighbors_file, index=False)
