#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 48:00:00
#SBATCH --job-name gather_ids
#SBATCH --partition RM
#SBATCH -o slurm-%x.out
#SBATCH --mem=200G


source ~/.bashrc
conda activate default

srun python3 gather_IDs.py 

