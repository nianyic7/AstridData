#!/bin/bash                                                                             
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=184,188,192,194,195,199,203,208,212,218,222,227,229,234,239,242,243,247,252,257,262,266,271,272,274,275,279 
#SBATCH --time 20:00:00
#SBATCH --job-name io-sort
#SBATCH --partition RM
#SBATCH -o slurm-%x-%a.out
#SBATCH --mem=100G

source ~/.bashrc
conda activate default

snap=$SLURM_ARRAY_TASK_ID
path="/hildafs/datasets/Asterix/BH_details_bigfile"

srun  python3 ../io_sorted_cols.py --snap $snap --path $path
