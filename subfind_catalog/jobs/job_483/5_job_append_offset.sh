#!/bin/bash 
#SBATCH -p development
#SBATCH --time=2:00:00
#SBATCH --nodes=8
#SBATCH -n 224
##SBATCH --mem=200GB
#SBATCH --job-name=wseed3
#SBATCH --output=slurm-append-offset.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nianyic@andrew.cmu.edu

module unload python3
source activate fast-mpi4py
which python

snap=483
dataroot="/home1/08942/nianyic/scratch3/Astrid"
subroot="$dataroot/subfind/subfind-$snap"
pigfile="/home1/08942/nianyic/ASTRID2_PIG/PIG_$snap"
dest="$dataroot/PIG_483_subfind"

mpirun -np 224 python3 ../code/4_append_suboff.py  --dest $dest --minpart 200 --subroot $subroot
