import re
import wget
import os
import gzip
import numpy as np
import pandas as pd
import subprocess
import csv
import glob
from collections import defaultdict,Counter
from models import Model


def create_csv_rasqual(path):

#     if os.path.exists(path+'all_rasqual.csv'):
#         print("rasqual dataset already read at {}\n".format(path+'all_rasqual.csv'))
#         return
    rasqual = gzip.open(path+'all.vars.inside.peak.gz')
    df_rasqual = pd.read_csv(rasqual,delimiter='\t',header=None)
    
    #indexes of columns that we need in the csv file
    peak,chromosome,SNP_pos,q_value,effect_size = 0,2,3,9,11 
    
    # add "chr" prefix to all entries in chromosomes column to match notation with other datasets
    chromosomes = list(df_rasqual.iloc[:,chromosome])
    for i in range(len(chromosomes)):
        chromosomes[i] = "chr" + str(chromosomes[i])
    df_rasqual.iloc[:,chromosome] = chromosomes
    
    # only keep columns that we need
    df_rasqual = df_rasqual.iloc[:,[peak,chromosome,SNP_pos,q_value,effect_size]]
    headers = ['peak','chr','start','Log_10 Benjamini-Hochberg Q-value','Effect size']
    df_rasqual.columns = headers
    df_rasqual['start'] = df_rasqual['start'] - 1
    df_rasqual['end'] =  df_rasqual['start'] + 1
    headers.append('end')
    df_rasqual.drop_duplicates(subset=['chr','start','end'],inplace=True)
    print(df_rasqual.shape)
    df_rasqual.to_csv(path+'all_rasqual.csv',sep='\t',header=headers,index=None)
    

def get_common_snps_compendium(path):
    # find all common snps bewteen the rasqual file and the compendium files. Include compendium.

#     if os.path.exists(path+'validation/rasqual/common_snp_compendium.csv'):
#         print("common snps already read at {}\n".format(path+'validation/rasqual/common_snp_compendium.csv'))
#         return

    df_rasqual = pd.read_csv(path+'validation/rasqual/all_rasqual.csv',delimiter='\t')
    df_rasqual.set_index(['chr','start','end'],inplace=True)

    df_compendium = pd.read_csv(path+'validation/all_snp_compendium.csv',delimiter='\t')
    df_compendium.set_index(['chr','start','end'],inplace=True)
    
    
    common_snps_compendium = df_compendium.merge(df_rasqual, left_index=True, right_index=True)
    common_snps_compendium.reset_index(inplace=True)

#     common_snps_index = df_compendium.index.intersection(df_rasqual.index)
    
#     common_snps_compendium = df_compendium.loc[common_snps_index,:]
#     common_snps_compendium['Log_10 Benjamini-Hochberg Q-value'] = df_rasqual['Log_10 Benjamini-Hochberg Q-value']
#     common_snps_compendium['Effect size'] = df_rasqual['Effect size']

    df_compendium.reset_index(inplace=True)
    df_rasqual.reset_index(inplace=True)
    print(common_snps_compendium.shape)
    print(common_snps_compendium.head())
    common_snps_compendium.to_csv(path+'validation/rasqual/common_snp_compendium.csv',index=False, sep='\t')


def get_common_snps_unique(path, FDR=10):
    # get all unique snps that are common bewteen the rasqual file and the compendium files. 
    # This only includes the chromosome and start position without the other columns
    
    file_path = path+'common_snp_unique_FDR={}%.csv'.format(FDR)
    if os.path.exists(file_path):
        print("common unique snps already exists at {}\n".format(file_path))
        return

    df_effectsize_labels = pd.DataFrame(columns=['Effect size'])

    df = pd.read_csv(path+'common_snp_compendium_FDR={}%.csv'.format(FDR), delimiter='\t')
    
    df.drop_duplicates(subset=['chr','start','end'], inplace=True)
    df.set_index(['chr','start','end'],inplace=True)
    df = df.sort_index()
    df.reset_index(inplace=True)
    df_valid_unique = df.iloc[:,0:3]
    df_effectsize_labels['Effect size'] = df['Effect size']
    print(df_valid_unique.shape)
    df_effectsize_labels.to_csv(path+'common_snp_effectsize_labels_FDR={}%.csv'.format(FDR), index=False, sep='\t')
    df_valid_unique.to_csv(path+'common_snp_unique_FDR={}%.csv'.format(FDR),index=False, sep='\t')
    


def filter_FDR(path, FDR=10):
    file_path = path+'common_snp_compendium_FDR={}%.csv'.format(FDR)
    if os.path.exists(file_path):
        print("common snps compendium filtered by FDR already exists at {}\n".format(file_path))
        return
    
    df = pd.read_csv(path+'common_snp_compendium.csv', delimiter='\t')
    # FDR of 10%
    df = df.loc[df['Log_10 Benjamini-Hochberg Q-value']<= np.log10(FDR/100.0)]
    print(df.shape)
    df.to_csv(path+'common_snp_compendium_FDR={}%.csv'.format(FDR),index=False, sep='\t')



def make_snp_footprint_matrix(path, flanking_size=0, FDR=10):
    # make snp vectors(A matrix) - a snp vector is a binary vector, each entry corresponds a specific motif.
    # Entry is 1 if that motifs that overlap with the snp, 0 otherwise. (Similar to how we built training matrix)
    motif_files = sorted(glob.glob(path+'combo/*.gz'))
    
    headers = ['chr','start','end']
    df = pd.read_csv(path+'validation/rasqual/common_snp_unique_FDR={}%.csv'.format(FDR),sep='\t')
    print("Making snp footprint matrix")
    
    '''
    if len(df.columns) > 3:
        print("snp footprint matrix already made at {}\n".format(path+'validation/dsQTL/common_snp_unique.csv'))
        return
    '''
        
    for idx,motif in enumerate(motif_files):
        motif_name = re.compile('.*/(.*).combo.bed.gz').search(motif).group(1)
        cmd_str = "bedtools intersect -a {}validation/rasqual/common_snp_unique_FDR={}%.csv -b {} -c | cut -f 4".format(path, FDR, motif)
        out = subprocess.check_output(cmd_str,shell=True)
        # convert numbers from string to integers. Last line is empty so don't include
        footprint = map(int,out.split('\n')[:-1])
        footprint = (np.array(footprint) >=1).astype(int)
        headers.append(motif_name)
        print("{}: {}, count of 1's: {}".format(motif_name, sum(footprint),Counter(footprint)[1]))
        df[motif_name] = footprint
        
        
    df.to_csv(path+'validation/rasqual/common_snp_unique_footprint_FDR={}%.csv'.format(FDR),index=False, header = headers, sep='\t')
    

# make reference/alternate matrix that is to be fed into the NN. Each row is a unique SNP with its motif footprint 
def make_ref_alt_matrix(path, allele='ref',flanking_size=0, FDR=10):
    print("Making {} allele matrix".format(allele))
    
#     if os.path.exists(path+'rasqual/common_snp_unique_{}.csv'.format(allele)):
#         print("{} matrix already made at {}\n".format(allele, path+'rasqual/common_snp_unique_{}.csv'.format(allele)))
#         return

    df_compendium = pd.read_csv(path+'rasqual/common_snp_compendium_FDR={}%.csv'.format(FDR), sep='\t')
    df_validation = pd.read_csv(path+'rasqual/common_snp_unique_footprint_FDR={}%.csv'.format(FDR),sep='\t')

    df_compendium.set_index(['chr','start','end'],inplace=True)
    df_validation.set_index(['chr','start','end'],inplace=True)
    
    # only interested in rows where the motif has an effect on binding (effect = 2)
    df_compendium = df_compendium.loc[df_compendium['effect']==2]

    ref_priorlodds = np.array(df_compendium['ref_priorlodds'])
    alt_priorlodds = np.array(df_compendium['alt_priorlodds'])
    if allele=='ref':
        # if ref_priorlodds-alt_priorlodds is greater than 0, reference allele increases binding
        df_compendium['binding_direction'] = (ref_priorlodds-alt_priorlodds > 0).astype(int)
    else:
        # if ref_priorlodds-alt_priorlodds less than 0, alternate allele increases binding
        df_compendium['binding_direction'] = (ref_priorlodds-alt_priorlodds <= 0).astype(int)
        
    for index in df_validation.index.values:
        try:
            centisnp = df_compendium.loc[index,['motif','binding_direction']]
        except:
            continue
        for motif in centisnp.values:
            motif_name, binding_direction = motif
            # only consider motifs that overlap with the snp
            if df_validation.at[index,motif_name] == 1:
                df_validation.at[index,motif_name] = binding_direction

    df_validation.reset_index(inplace=True)
    df_compendium.reset_index(inplace=True)
    df_validation.to_csv(path+'rasqual/common_snp_unique_{}_FDR={}%.csv'.format(allele,FDR),index=False, sep='\t')



def log_odds_diff(p1,p2):
    np.seterr(all='raise') 
    # avoid taking the log of 0
    p1[p1==1.0] = 0.999999
    p2[p2==1.0] = 0.999999
    p1[p1==0.0] = 0.000001
    p2[p2==0.0] = 0.000001
    return np.log(p1) - np.log(1-p1) - (np.log(p2) - np.log(1-p2))

def make_circuitSNP_predictions(path, seed=0, FDR=10):
    INPUT_DIM = 1372
    nn_models = Model(INPUT_DIM)
    
    dir_path = '{}data/validation/rasqual/results/FDR={}%/models_seed={}/'.format(path,FDR,seed)
    
    # make directory to save results to 
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)
    
    models_dir= sorted(glob.glob('{}models/models_seed={}/*.h5'.format(path,seed)))
    
    X_data_ref = pd.read_csv(path+'data/validation/rasqual/common_snp_unique_ref_FDR={}%.csv'.format(FDR), delimiter='\t')
    X_data_alt = pd.read_csv(path+'data/validation/rasqual/common_snp_unique_alt_FDR={}%.csv'.format(FDR), delimiter='\t')

    output = pd.DataFrame(columns = ['ref','alt','circuitsnp_prediction'])
                             
    for model_weights in models_dir:
        print(model_weights)
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
        output.to_csv(dir_path+'model{}_replicate_{}_validation_output.csv'
                      .format(model_arch, replicate_num),index=False,sep='\t')
        

def driver(ROOT_DIR):
    path = ROOT_DIR+'data/validation/rasqual/'
    if not os.path.exists(path+'all.vars.inside.peak.gz'):
        wget.download('http://genome.grid.wayne.edu/centisnps/rasqual/all.vars.inside.peak.gz',path+'all.vars.inside.peak.gz')
    
    # create_csv_rasqual(path) #519013
    # get_common_snps_compendium(ROOT_DIR+'data/')  #638169


    results_path = path + 'results/'
    if not os.path.exists(results_path):
        os.mkdir(results_path)

    for fdr in (10,5):
        filter_FDR(path, FDR=fdr)

        get_common_snps_unique(path, FDR=fdr)

        make_snp_footprint_matrix(ROOT_DIR+'data/',flanking_size=0, FDR=fdr)
        make_ref_alt_matrix(ROOT_DIR+'data/validation/',allele='ref', flanking_size=0, FDR=fdr)
        make_ref_alt_matrix(ROOT_DIR+'data/validation/',allele='alt', flanking_size=0, FDR=fdr)  

                
        if not os.path.exists(results_path+'FDR={}%'.format(fdr)):
            os.mkdir(results_path+'FDR={}%'.format(fdr))

        for i in range(0,10):
            make_circuitSNP_predictions(ROOT_DIR, seed=i, FDR=fdr)