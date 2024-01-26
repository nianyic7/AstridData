#!/bin/bash 
#SBATCH -p RM
#SBATCH --time=00:20:00
#SBATCH --nodes=2
#SBATCH --mem=100GB
#SBATCH --job-name=wseed3
#SBATCH --output=slurm-test-mpi.out

source activate default

mpirun -n 56 python3 fileMPI_test.py
