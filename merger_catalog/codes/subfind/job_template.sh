#!/bin/bash 
#SBATCH -p RM
#SBATCH --time=5:00:00
#SBATCH --nodes=1
##SBATCH -n 128
##SBATCH --ntasks-per-node=8
#SBATCH --mem=100GB
#SBATCH --job-name=sub-snap
#SBATCH --output=slurm-sub-snap.out

snap=sub-snap
WORKDIR="/hildafs/home/nianyic/asterix_bh/subgroups"

srun -n 1 python3 groups2hdf5_z34.py --snap $snap &> $WORKDIR/sub-snap-out.out

