#!/bin/bash
#SBATCH -p RM
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --job-name=extract
#SBATCH --time=20:00:00
#SBATCH --mem=200GB
#SBATCH --output=slurm-group2hdf5.log

savedir="/hildafs/home/nianyic/scratch1/Astrid_data/binary_cat/subfind/hdf5Groups/"
for snap in 272 282
do
srun python3 groups2hdf5.py --snap $snap  --savedir $savedir 
done
