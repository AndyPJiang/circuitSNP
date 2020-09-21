import numpy as np
import glob
import os
import pandas as pd
from models import Model
import re


# Calculate the log odds difference
def log_odds_diff(p1,p2):
    np.seterr(all='raise') 
    # avoid taking the log of 0
    p1[p1==1.0] = 0.999999
    p2[p2==1.0] = 0.999999
    p1[p1==0.0] = 0.000001
    p2[p2==0.0] = 0.000001
    
    return np.abs(np.log(p1) - np.log(1-p1) - (np.log(p2) - np.log(1-p2)))


# generate circuitSNP prediction by taking the log odds difference between the reference matrix
# output and alternate matrix output
def make_circuitSNP_predictions(path,seed,flanking_size):
    INPUT_DIM = 1372
    nn_models = Model(INPUT_DIM)
    
    # make directory to save results to 
    if not os.path.exists('{}data/validation/dsQTL/results/flanking_size={}/models_seed={}'.format(path, flanking_size, seed)):
        os.mkdir('{}data/validation/dsQTL/results/flanking_size={}/models_seed={}'.format(path, flanking_size, seed))
    
    models_dir= sorted(glob.glob('{}models/models_seed={}/*.h5'.format(path,seed)))
    
    X_data_ref = pd.read_csv(path+'data/validation/dsQTL/validation_matrices_flanksize={}/common_snp_unique_ref.csv'.format(flanking_size), delimiter='\t')
    X_data_alt = pd.read_csv(path+'data/validation/dsQTL/validation_matrices_flanksize={}/common_snp_unique_alt.csv'.format(flanking_size), delimiter='\t')

    output = pd.DataFrame(columns = ['ref','alt','circuitsnp_prediction'])
                             
    for model_weights in models_dir:
        model_arch = re.compile('.*/model(.*)_replicate.*').search(model_weights).group(1)
        model_arch_split = model_arch.split('+')
        replicate_num = re.compile('/.*_replicate_(.*).h5').search(model_weights).group(1)

        hid_layer1_units, hid_layer2_units = map(int,(model_arch_split[0][1:],model_arch_split[1][:-1]))

        if hid_layer1_units == 0 and hid_layer2_units == 0:
            model = nn_models.create_model_simple0()
        elif hid_layer2_units == 0:
            model = nn_models.create_model_simple1()
        else:
            model = nn_models.create_model(hid_layer1_units,hid_layer2_units)

        model.load_weights(model_weights)
        output['ref'] = model.predict(X_data_ref.iloc[:,3:],verbose=0).flatten()
        output['alt'] = model.predict(X_data_alt.iloc[:,3:],verbose=0).flatten()
        
        
        output['circuitsnp_prediction'] = log_odds_diff(np.array(output['ref']),np.array(output['alt']))
        output.to_csv('{}data/validation/dsQTL/results/flanking_size={}/models_seed={}/model{}_replicate_{}_validation_output.csv'
                      .format(path, flanking_size, seed, model_arch,replicate_num),index=False,sep='\t')
        

def driver(ROOT_DIR,seed=42,flanking_size=0):
    if not os.path.exists('{}data/validation/dsQTL/results/flanking_size={}'.format(ROOT_DIR, flanking_size)):
        os.mkdir('{}data/validation/dsQTL/results/flanking_size={}'.format(ROOT_DIR, flanking_size))
        
    print("Making predictions for models with seed = {}".format(seed))
    make_circuitSNP_predictions(ROOT_DIR,seed, flanking_size)
