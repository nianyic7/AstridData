#!/bin/bash                                                                             
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 40:00:00
#SBATCH --job-name sort-734-751
#SBATCH --partition HENON
#SBATCH -o slurm-%x.out
#SBATCH --mem=400G

source ~/.bashrc
conda activate default

sbeg=734
send=751
path="/hildafs/datasets/Asterix/BH_details_bigfile2"

srun  python3 ../sort_and_count.py --sfirst $sbeg --slast $send --path $path
