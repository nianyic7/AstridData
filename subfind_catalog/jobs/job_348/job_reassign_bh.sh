#!/bin/bash 
#SBATCH -p RM
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH --mem=100GB
#SBATCH --job-name=wseed3
#SBATCH --output=slurm-reassign-bh.out


srun python3 ../reassign_bh.py  --snap 348
