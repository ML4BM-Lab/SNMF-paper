
from pathlib import Path

import torch
torch.set_num_threads(16)

import sys
sys.path.append("/scratch/lalonsoeste/PhD/NMF_deconvolution/methods/SpiceMix/repo/SpiceMix")
from model import SpiceMix

import pandas as pd
import numpy as np

from sklearn.neighbors import kneighbors_graph

import pickle


context = dict(device='cuda' if torch.cuda.is_available() else 'cpu', dtype=torch.float64)

data_path = Path(sys.argv[1])
output_path = Path(sys.argv[2])
K = int(sys.argv[3])

result_path = output_path / 'tmp' / 'SpiceMix.h5'

obj = SpiceMix(
    K=K,
    lambda_Sigma_x_inv=1e-6, 
    power_Sigma_x_inv=2,
    repli_list=['r1'],
    context=context,
    context_Y=context,
    path2result=result_path,
)

def load_from_csv(self, csv_path, replicate_name='r1'):
    expr = pd.read_csv(csv_path, index_col=0).T

    self.repli_list = [replicate_name]  # mock one replicate
    self.Ys = [torch.tensor(expr.values, **self.context_Y)]
    self.Ns, self.Gs = [expr.shape[0]], [expr.shape[1]]
    self.GG = self.Gs[0]
    self.genes = [list(expr.columns)]
    
    spot_locations = pd.DataFrame(expr.index.str.split("x").tolist(), columns=["array_row", "array_col"])
    spot_locations = spot_locations[["array_row", "array_col"]].values.astype(np.float32)

    E = np.array(kneighbors_graph(spot_locations, n_neighbors=20, mode='connectivity', include_self=False).todense())
    
    self.Es = [E]
    self.Es_isempty = [sum(map(len, E)) == 0]
    
    expr['repli'] = replicate_name
    self.meta = expr

load_from_csv(obj, data_path)

# obj.load_dataset = load_from_csv.__get__(obj)
# obj.load_dataset(data_path)

with open(output_path / "tmp" / "spicemix_obj.pickle", 'wb') as handle:
    pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)