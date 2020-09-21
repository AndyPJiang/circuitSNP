#!/bin/bash
#SBATCH --nodes 6
#SBATCH --account=punim0614
#SBATCH --qos=gpgpumdhs
#SBATCH --partition gpgpu
#SBATCH --gres=gpu:p100:4
#SBATCH --time 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --job-name=train_models
#SBATCH --mem=100000

#load modules
module load Python/2.7.13-GCC-6.2.0
module load Tkinter/2.7.11-intel-2017.u2-Python-2.7.11
module load Python/2.7.13-GCC-6.2.0
module load Tensorflow/1.11.0-intel-2017.u2-Python-2.7.13-GPU
module load SAMtools/1.8-GCC-6.2.0-HTSlib-1.8
module load BEDTools/2.27.1-spartan_gcc-6.2.0
module load web_proxy
source ${HOME}/p27_jp/bin/activate

#execute relevant python script
python train.py
