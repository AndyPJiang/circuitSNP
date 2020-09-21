from bs4 import BeautifulSoup
from requests import get
import re
import wget
import os
import urllib
import gzip
import tarfile
import csv
import glob


def get_combo_files_list(path, URL_COMBO, TOTAL_EXPECTED_COMBO_FILES):
    #Make a text file containing the TF motif datasets. The datasets are used to construct the input to the NN
    try:
        response = get(URL_COMBO)
        soup = BeautifulSoup(response.text, 'html.parser')

        files = soup.find_all('a', string = re.compile('.*.combo.bed.gz.*'))
        f = open(path+'combo_files.txt','w')

        count = 0
        for file in files:
            count+=1
            f.write(file.string + '\n')

        print("Total files found: "+ str(count))
        print("Expected files: " + str(TOTAL_EXPECTED_COMBO_FILES))
        f.close()
    except:
        print("Unable to get files")

def get_compendium_files_list(path, URL_COMPENDIUM, TOTAL_EXPECTED_COMPENDIUM_FILES):
    # Make a text file containing the SNP locations and intersections
    try:
        response = get(URL_COMPENDIUM)
        soup = BeautifulSoup(response.text, 'html.parser')

        files = soup.find_all('a', string = re.compile('.*.bed.gz'))
        f = open(path+'compendium_files.txt','w')

        count = 0
        for file in files:
          count+=1
          f.write(file.string + '\n')

        print("Total files found: "+ str(count))
        print("Expected files: " + str(TOTAL_EXPECTED_COMPENDIUM_FILES))

        f.close()
    except:
        print("Unable to get files")


def make_all_directories(ROOT_DIR):
    # create data directories if they don't exist already
    if not os.path.exists(ROOT_DIR+'data/'):
        os.mkdir(ROOT_DIR+'data/')
    if not os.path.exists(ROOT_DIR+'data/combo/'):
        os.mkdir(ROOT_DIR+'data/combo')
    if not os.path.exists(ROOT_DIR+'data/compendium/'):
        os.mkdir(ROOT_DIR+'data/compendium/')
    if not os.path.exists(ROOT_DIR+'data/validation/'):
        os.mkdir(ROOT_DIR+'data/validation/')
    if not os.path.exists(ROOT_DIR+'models/'):
        os.mkdir(ROOT_DIR+'models/')
    if not os.path.exists(ROOT_DIR+'data/validation/results/'):
        os.mkdir(ROOT_DIR+'data/validation/results/')
    if not os.path.exists(ROOT_DIR+'data/validation/dsQTL/'):
        os.mkdir(ROOT_DIR+'data/validation/dsQTL/')
    if not os.path.exists(ROOT_DIR+'data/validation/rasqual/'):
        os.mkdir(ROOT_DIR+'data/validation/rasqual/')

#function used to download TF motif and SNP dataset
def download_dataset(path, to_download, downloaded, failed_download,url):
    # all files to download
    f_to_download = open(to_download,'r')

    # files already downloaded. 
    try:
        f_downloaded = open(downloaded,'r')
    except:
        f_downloaded = open(downloaded,'w+')

    # files that failed to download
    f_failed_download = open(failed_download,'a+')

    # get all the files that already have been downloaded as a set for constant time lookup
    already_downloaded = set(map(str.strip,f_downloaded.readlines()))

    f_downloaded.close()

    f_downloaded = open(downloaded,'a+')

    success = 0
    failed = 0
    count_tries = 0

    for file in f_to_download.readlines():
        file = file.strip()
        count_tries +=1
        # don't download files that already exist
        if file in already_downloaded:
            #print("File " + str(count_tries) + ' already downloaded')
            continue
        try:
            wget.download(url + file, out=path)
            f_downloaded.write(file + '\n')
            #print("Successfully downloaded file " + str(count_tries) +' '+ file)
            success+=1
        except:
            f_failed_download.write(file + '\n')
            #print("Failed to download file " + str(count_tries) +' '+ file)
            failed +=1

    print("Completed")
    print("Summary:")
    print("Total number of files attempted: " + str(count_tries))
    print(str(success)+' files downloaded successfully')
    print(str(failed)+' files failed to download\n')

    f_to_download.close()
    f_downloaded.close()
    f_failed_download.close()
    
    
#function used to download beer labs chromatin windows dataset
def download_windows_dataset(path, url):
    wget.download(url, out=path)
    in_file = tarfile.open(path+'gm12878_sequence_sets.tar.gz', 'r')
    in_file.extractall(path)
    
    # remove gzip file after downloading and unzipping
    os.remove(path+'gm12878_sequence_sets.tar.gz')


def extract_windows(path, files,expected_windows_found, type='open'):
    windows = []
    
    for file in files:
        f = Fasta(file)
        for window in sorted(f.keys()):
            chrsm_num, window_range = window.split(':')
            [start, end] = map(int,window_range.split('-'))
            windows.append([chrsm_num,start-1,end])
    #windows.sort(key = lambda x:(int(x[0][3:]),x[1],x[2]))
    if os.path.exists(path+'dnase_windows_{}.bed'.format(type)):
        print("File Already exists. Overwriting...")

    df = pd.DataFrame(windows)
    df.to_csv(path + 'dnase_windows_{}.bed'.format(type), index=False, header=None,sep='\t')
    
    print("Total of {} {} dnase windows found".format(len(windows),type))
    print("Expected {} {} dnase windows found\n".format(expected_windows_found,type))




def download_dnase_dataset(path, URL_WINDOWS,TOTAL_EXPECTED_OPEN_WINDWOS, TOTAL_EXPECTED_CLOSED_WINDWOS):
    print("Downloading Dnase Windows Dataset...")
    if os.path.exists(path+'input_windows.bed'):
        print("Dnase Windows Dataset already downloaded, aborting...\n")
        return
        
    # download beers lab dataset for dnase windows, extract data and concatenate into one file
    download_windows_dataset(path, URL_WINDOWS)
    
    window_data_dir = path+'gm12878_sequence_sets/'
    
    # all files that contain data on open windows
    open_files=[window_data_dir+'gm12878_shared.fa']
    # all files that contain data on closed windows. Ratio of dataset is 1 open : 5 closed 
    closed_files=[window_data_dir + 'nullseqs_gm12878_shared.{}.1.fa'.format(i) for i in range(1,6)]

    # extract the windows separately into csv files
    extract_windows(path, open_files, TOTAL_EXPECTED_OPEN_WINDWOS,type='open')
    extract_windows(path, closed_files,TOTAL_EXPECTED_CLOSED_WINDWOS,type='closed')

    print("Combining closed and open windows into single file...")
    # contatentate both open and closed windows into one file
    df_open = pd.read_csv(path+'dnase_windows_open.bed',sep='\t',header=None)
    df_closed = pd.read_csv(path+'dnase_windows_closed.bed',sep='\t',header=None)
    frame = pd.concat([df_open,df_closed], axis=0, ignore_index=True)
    frame.to_csv(path+'input_windows.bed', index=False,header=None,sep='\t')
    
    
    print("Download Complete\n")
    # print first five rows of contatenated file
    #df = pd.read_csv(ROOT_DIR+'data/input_windows.bed')
    #print("First 5 rows of concatenated csv file of dnase windows")
    #df.head()



def driver(ROOT_DIR):
    URL_COMBO = 'http://genome.grid.wayne.edu/centisnps/combo/'
    TOTAL_EXPECTED_COMBO_FILES = 1633

    URL_COMPENDIUM = 'http://genome.grid.wayne.edu/centisnps/compendium/'
    TOTAL_EXPECTED_COMPENDIUM_FILES = 1371

    URL_WINDOWS = 'http://www.beerlab.org/deltasvm/downloads/gm12878_sequence_sets.tar.gz'
    TOTAL_EXPECTED_OPEN_WINDWOS = 22384
    TOTAL_EXPECTED_CLOSED_WINDWOS = 111920

    

    #get_combo_files_list(ROOT_DIR,URL_COMBO,TOTAL_EXPECTED_COMBO_FILES)
    #get_compendium_files_list(ROOT_DIR,URL_COMPENDIUM,TOTAL_EXPECTED_COMPENDIUM_FILES)
    make_all_directories(ROOT_DIR)
    
    # download TF motif(combo) dataset
    path_combo = ROOT_DIR+'/data/combo/'
    all_files = 'combo_files.txt'
    already_downloaded = 'combo_downloaded.txt'
    failed_download = 'combo_failed_download.txt'
    print("Downloading Combo Dataset...")
    download_dataset(path_combo, all_files, already_downloaded, failed_download, URL_COMBO)

    # download SNP(compendium) dataset
    path_compendium = ROOT_DIR+'/data/compendium/'
    all_files = 'compendium_files.txt'
    already_downloaded = 'compendium_downloaded.txt'
    failed_download = 'compendium_failed_download.txt'
    print("Downloading Compendium Dataset...")
    download_dataset(path_compendium, all_files, already_downloaded, failed_download, URL_COMPENDIUM)
    
    # Download dnase windows dataset
    download_dnase_dataset(ROOT_DIR+'data/',URL_WINDOWS,TOTAL_EXPECTED_OPEN_WINDWOS, TOTAL_EXPECTED_CLOSED_WINDWOS)
    
    print("Download Finished")