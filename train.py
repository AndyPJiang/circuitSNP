# slurm job id: 19465764

import pandas as pd
import numpy as np
from modules import train_models
from modules import make_circuitSNP_predictions


ROOT_DIR = "/data/projects/punim0614/andy/circuitSNP/"
df = pd.read_csv(ROOT_DIR+'data/train_data_mat.csv',sep='\t')


for seed in range(10):
    train_models.driver(ROOT_DIR, df, seed=seed)
    #make_circuitSNP_predictions.driver(ROOT_DIR, seed)
