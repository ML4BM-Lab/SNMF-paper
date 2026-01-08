
from pathlib import Path

import sys

import torch
torch.set_num_threads(16)

import numpy as np
import pandas as pd

import pickle

output_path = Path(sys.argv[1])
niter = int(sys.argv[2])
seed = int(sys.argv[3])

sys.path.append("/scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SpiceMix/repo/SpiceMix")
with open(output_path / "tmp" / "spicemix_obj.pickle", 'rb') as handle:
    obj = pickle.load(handle)

result_path = output_path / 'tmp' / 'SpiceMix.h5'

np.random.seed(seed)

num_pcs = 50
n_neighbors = 20
res_lo = .1
res_hi = 5.

obj.initialize(
    method='louvain', kwargs=dict(num_pcs=num_pcs, n_neighbors=n_neighbors, resolution_boundaries=(res_lo, res_hi), num_rs=30),
)

for iiter in range(10):
    obj.estimate_weights(iiter=iiter, use_spatial=[False])
    obj.estimate_parameters(iiter=iiter, use_spatial=[False])
obj.initialize_Sigma_x_inv()
for iiter in range(1, niter+1):
    print(f'Iteration {iiter}')
    obj.estimate_parameters(iiter=iiter, use_spatial=[True])
    obj.estimate_weights(iiter=iiter, use_spatial=[True])

import h5py
with h5py.File(result_path, 'r') as f:
    proportions = pd.DataFrame(f[f"latent_states/XT/r1/{niter}"][()], index=obj.meta.index)
    metagenes = pd.DataFrame(f[f"parameters/M/{niter}"][()], index=obj.genes)

proportions = proportions.div(proportions.sum(axis=1), axis=0)
proportions.to_csv(output_path / "SpiceMix_proportions.csv")
metagenes.to_csv(output_path / "SpiceMix_metagenes.csv")
