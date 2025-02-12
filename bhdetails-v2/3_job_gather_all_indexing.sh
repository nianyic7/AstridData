#!/bin/bash                                                                             
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 10:00:00
#SBATCH --job-name save-idx
#SBATCH --partition RM
#SBATCH -o slurm-%x.out
#SBATCH --mem=200G

source ~/.bashrc
conda activate default


srun  python3 3_gather_all_indexing.py
