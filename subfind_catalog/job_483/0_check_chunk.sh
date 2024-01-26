#!/bin/bash 
#SBATCH -p small
#SBATCH --time=3:00:00
#SBATCH --nodes=2
##SBATCH --mem=200GB
#SBATCH --ntasks-per-node=14
#SBATCH --job-name=save-order
#SBATCH --output=slurm-check-chunk.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nianyic@andrew.cmu.edu

module unload python3
source activate fast-mpi4py 
which python
snap=483
dataroot="/home1/08942/nianyic/scratch3/Astrid"
tabfile="fof_subhalo_tab_$snap.hdf5"
grpfile="snap_$snap.hdf5"
subroot="$dataroot/subfind/subfind-$snap"
savedir="$dataroot/subfind/checks-$snap"

ibrun -n 28 python3 ../code/0_check_chunk.py --subroot $subroot --tabfile $tabfile --grpfile $grpfile --savedir $savedir 
