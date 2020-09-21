# slurm job id: 19505258
# make validation predictions 

import os
from modules import make_circuitSNP_predictions_dsQTL
ROOT_DIR = "/data/projects/punim0614/andy/circuitSNP/"

flank_sizes = (0,10,100,300)

if not os.path.exists('{}data/validation/dsQTL/results/'.format(ROOT_DIR)):
    os.mkdir('{}data/validation/dsQTL/results/'.format(ROOT_DIR))
        
for flank_size in flank_sizes:
    for seed in range(10):
        make_circuitSNP_predictions_dsQTL.driver(ROOT_DIR, seed=seed, flanking_size=flank_size)
