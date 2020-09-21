import pandas as pd
import subprocess
import numpy as np
import os
import glob
import re
from collections import Counter

def make_training_data(path):
    motif_files = sorted(glob.glob(path+'combo/*.gz'))
    headers = ['chrsm','start','end']
    df = pd.read_csv(path+'input_windows.bed',sep='\t',names = ['chrsm','start','end'])
    print("Making training data...")
    if os.path.exists(path+'train_data_mat.csv'):
        print("Training data already made at {}".format(path+'train_data_mat.csv'))
        return

    for idx,motif in enumerate(motif_files):
        motif_name = re.compile('.*/(.*).combo.bed.gz').search(motif).group(1)
        cmd_str = "bedtools intersect -a data/input_windows.bed -b {} -c | cut -f 4".format(motif)
        out = subprocess.check_output(cmd_str,shell=True)
        # convert numbers from string to integers. Last line is empty so don't include
        footprint = map(int,out.split('\n')[:-1])
        footprint = (np.array(footprint) >=1).astype(int)
        headers.append(motif_name)
        print("{}: {}, count of 1's: {}".format(motif_name, sum(footprint),Counter(footprint)[1]))
        df[motif_name] = footprint
    df.to_csv(path+'train_data_mat.csv',index=False, header = headers, sep='\t')
    

def make_training_labels(path):
    df = pd.read_csv(path+'train_data_mat.csv',sep='\t')
    print("Making training labels...")
    if 'label' in df.columns:
        print("Labels already made at {}".format(path+'train_data_mat.csv'))
        return
    df_open = pd.read_csv(path+'dnase_windows_open.bed', header=None)
    df_closed = pd.read_csv(path+'dnase_windows_closed.bed', header=None)
    
    labels_open = np.ones((df_open.size,), dtype=int)
    labels_closed = np.zeros((df_closed.size,), dtype=int)
    all_labels = np.concatenate([labels_open,labels_closed], axis=None)
    df['label'] = all_labels
    df.to_csv(path+'train_data_mat.csv',index=False, sep='\t')


def driver(ROOT_DIR):
    make_training_data(ROOT_DIR+'data/')
    make_training_labels(ROOT_DIR+'data/')
