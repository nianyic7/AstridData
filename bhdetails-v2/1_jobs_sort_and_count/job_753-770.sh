#!/bin/bash                                                                             
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 40:00:00
#SBATCH --job-name sort-753-770
#SBATCH --partition HENON
#SBATCH -o slurm-%x.out
#SBATCH --mem=400G

source ~/.bashrc
conda activate default

sbeg=753
send=770
path="/hildafs/datasets/Asterix/BH_details_bigfile2"

srun  python3 ../sort_and_count.py --sfirst $sbeg --slast $send --path $path
