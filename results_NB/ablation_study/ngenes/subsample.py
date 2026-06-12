
import sys

import pandas as pd
import numpy as np

if __name__ == "__main__":
    data_path = sys.argv[1]
    markers_path = sys.argv[2]
    ngenes = int(sys.argv[3])
    output_path = sys.argv[4]
    mgs_output_path = sys.argv[5]

    data = pd.read_csv(data_path, index_col=0).T
    markers = pd.read_csv(markers_path)

    gene_var = data.var(axis=0)
    top_genes = gene_var.nlargest(ngenes).index

    downsampled = data[top_genes]
    downsampled.T.to_csv(output_path)

    new_markers = markers[markers["gene"].isin(top_genes)]

    if set(markers["cluster"].unique()) != set(new_markers["cluster"].unique()):
        all_genes = markers["gene"].unique().tolist()
        
        for cluster in set(markers["cluster"].unique()).difference(set(new_markers["cluster"].unique())):
            new_markers = pd.concat(
                [
                    new_markers, 
                    pd.DataFrame({
                        "cluster": [cluster] * len(all_genes), 
                        "gene": all_genes
                    })
                ], 
                axis=0, 
                ignore_index=True
            )

    new_markers.to_csv(mgs_output_path)


