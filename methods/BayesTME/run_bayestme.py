
import sys
import pickle

from bayestme.svi.deconvolution import deconvolve

import numpy as np

import pandas as pd
import os

output_path = sys.argv[1]
K = int(sys.argv[2])
rho = float(sys.argv[3])
seed = int(sys.argv[4])

with open(os.path.join(output_path, "tmp", "stdata.pickle"), 'rb') as handle:
    stdata = pickle.load(handle)
with open(os.path.join(output_path, "tmp", "obsnames.pickle"), 'rb') as handle:
    obs_names = pickle.load(handle)

rng = np.random.default_rng(seed=seed)
deconvolution_result = deconvolve(stdata=stdata, 
                                  n_components=K,
                                  rho=rho,
                                  rng=rng)

proportions = pd.DataFrame(deconvolution_result.cell_prob_trace.mean(axis=0), index=obs_names)
proportions.to_csv(os.path.join(output_path, "BayesTME_proportions.csv"))
