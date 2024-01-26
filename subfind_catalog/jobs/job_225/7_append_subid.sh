#!/bin/bash 
#SBATCH -p development
#SBATCH --time=2:00:00
#SBATCH --nodes=2
#SBATCH -n 56
##SBATCH --mem=100GB
#SBATCH --job-name=sidx
#SBATCH --output=slurm-app-sidx.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nianyic@andrew.cmu.edu

module unload python3
source activate fast-mpi4py
which python

snap=483
dataroot="/home1/08942/nianyic/scratch3/Astrid"
pigfile="/home1/08942/nianyic/ASTRID2_PIG/PIG_$snap"
outdir="$dataroot/PIG_483_subfind"
ibrun -n 56 python3 ../code/6_bh-subid.py  --dest $outdir --gstart 0


