#!/bin/bash 
#SBATCH -p RM
#SBATCH --time=2:00:00
#SBATCH --nodes=2
#SBATCH --mem=100GB
#SBATCH --job-name=wseed
#SBATCH --output=slurm-fileMPI.out



mpirun -n 56 python3 fileMPI_test.py
