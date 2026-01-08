
import sys
import os

import pickle

import pandas as pd
import scanpy as sc

import torch

from starfysh import utils
from starfysh import starfysh as sf_model

output_path = sys.argv[1]
learning_rate = float(sys.argv[2])
seed = int(sys.argv[3])

with open(os.path.join(output_path, "tmp", "visium_args.pickle"), 'rb') as handle:
    visium_args = pickle.load(handle)
adata = sc.read_h5ad(os.path.join(output_path, "tmp", "adata.h5ad"))
gene_sig = pd.read_csv(os.path.join(output_path, "tmp", "gene_sig.csv"), index_col=0)

## Run Starfysh
n_repeats = 4 # recommend > 3 for selecting a better trained model 
epochs = 200
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

model, loss = utils.run_starfysh(visium_args,
                                 n_repeats=n_repeats,
                                 epochs=epochs,
                                 lr=learning_rate,
                                 device=device,
                                 seed=seed)

inference_outputs, generative_outputs = sf_model.model_eval(model,
                                                            visium_args.adata_norm,
                                                            visium_args,
                                                            device=device)

prop_pred_df = pd.DataFrame(inference_outputs['qc_m'].cpu(), index=visium_args.adata_norm.obs_names, columns=gene_sig.columns)
prop_pred_df.to_csv(os.path.join(output_path, "starfysh_proportions.csv"))

