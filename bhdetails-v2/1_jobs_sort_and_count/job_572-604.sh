#!/bin/bash                                                                             
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 20:00:00
#SBATCH --job-name sort-572-604
#SBATCH --partition RM
#SBATCH -o slurm-%x.out
#SBATCH --mem=200G

source ~/.bashrc
conda activate default

sbeg=572
send=604
path="/hildafs/datasets/Asterix/BH_details_bigfile2"

srun  python3 ../1_sort_and_count.py --sfirst $sbeg --slast $send --path $path