
import pandas as pd
import numpy as np

from pathlib import Path

if __name__ == "__main__":
    data_dir = Path(__file__).resolve().parents[1]
    zenodo_dir = Path(__file__).resolve().parents[3] / "data" / "TNBC" / "final"
    output_dir = data_dir / "subsets"
    output_dir.mkdir(exist_ok=True)

    counts = pd.read_csv(zenodo_dir / "TNBC_counts.csv", index_col=0)
    proportions = pd.read_csv(zenodo_dir / "TNBC_proportions.csv", index_col=0)

    locations = pd.DataFrame([[int(v) for v in col] for col in counts.columns.str.split("x")], columns=["x","y"], index=counts.columns)
    center = locations.mean().values

    for N in [100, 500, 1000, 2000]:
        side = np.sqrt(N)
        n = 0
        while n < N:
            x_mask = (locations['x'] > center[0] - side/2) & (locations['x'] < center[0] + side/2)
            y_mask = (locations['y'] > center[1] - side/2) & (locations['y'] < center[1] + side/2)
            subset_spots = locations[x_mask & y_mask].index
            n = len(subset_spots)
            side += 2
        
        folder_name = output_dir / f"{n}spots_{side:.1f}x{side:.1f}"
        folder_name.mkdir(exist_ok=True)

        subset_counts = counts.loc[:,subset_spots]
        subset_proportions = proportions.loc[subset_spots,:]

        subset_counts.to_csv(folder_name / "counts.csv")
        subset_proportions.to_csv(folder_name / "proportions.csv")


