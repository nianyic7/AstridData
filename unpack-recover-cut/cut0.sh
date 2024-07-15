#!/bin/bash                                                                             
#SBATCH -p development
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=cut4-0
#SBATCH --output=%x-%j.out
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nianyic@andrew.cmu.edu

module unload python3
source activate fast-mpi4py
 
hostname; pwd; date
 
#export OMP_NUM_THREADS=1
 
ibrun -n 80  python gascut_mpi.py || exit 1
 
date
