#!/bin/bash
#SBATCH -p RM
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --job-name=extract
#SBATCH --time=20:00:00
#SBATCH --mem=200GB
#SBATCH --output=slurm-grouplists.log

savedir="/hildafs/home/nianyic/scratch1/Astrid_data/binary_cat/subfind/grpList"
for snap in 272 282
do
srun python3 merger_groups.py --snap $snap  --savedir $savedir 
done
