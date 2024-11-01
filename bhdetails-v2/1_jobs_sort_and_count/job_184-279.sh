#!/bin/bash                                                                             
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 8:00:00
#SBATCH --job-name sort-184-279
#SBATCH --partition RM
#SBATCH -o slurm-%x.out
#SBATCH --mem=200G

source ~/.bashrc
conda activate default

sbeg=184
send=279
path="/hildafs/datasets/Asterix/BH_details_bigfile"

srun  python3 ../1_sort_and_count.py --sfirst $sbeg --slast $send --path $path