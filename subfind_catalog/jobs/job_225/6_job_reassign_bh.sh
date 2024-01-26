#!/bin/bash 
#SBATCH -p small
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH -n 1
##SBATCH --mem=100GB
#SBATCH --job-name=reassifn
#SBATCH --output=slurm-reassign-bh.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nianyic@andrew.cmu.edu

module unload python3
source activate fast-mpi4py
which python

snap=483
dataroot="/home1/08942/nianyic/scratch3/Astrid"
pigfile="/home1/08942/nianyic/ASTRID2_PIG/PIG_$snap"
indir="$dataroot/PIG_483_subfind"
outdir="$dataroot/PIG_483_subfind_reassign"
srun python3 ../code/5_reassign_bh.py  --indir $indir --outdir $outdir --min_bhmass1 1e8 --min_gal_bh_ratio1 30 --min_bhmass2 1e7 --min_gal_bh_ratio2 2.4


