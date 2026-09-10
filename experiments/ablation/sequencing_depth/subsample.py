
import sys

import pandas as pd
import numpy as np

if __name__ == "__main__":
    data_path = sys.argv[1]
    depth_fraction = float(sys.argv[2])
    output_path = sys.argv[3]

    data = pd.read_csv(data_path, index_col=0).T

    rng = np.random.default_rng(42)

    downsampled = {}

    for gene in data.columns:
        counts = data[gene].values
        total_counts = counts.sum()

        if total_counts == 0:
            continue

        target_depth = int(total_counts * depth_fraction)

        if target_depth == 0:
            downsampled[gene] = np.zeros_like(counts)
            continue

        probs = counts / total_counts

        new_counts = rng.multinomial(target_depth, probs)

        downsampled[gene] = new_counts

    pd.DataFrame(
        downsampled,
        index=data.index
    ).T.to_csv(output_path)

