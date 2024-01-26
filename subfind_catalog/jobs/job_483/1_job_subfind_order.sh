#!/bin/bash 
#SBATCH -p development
#SBATCH --time=2:00:00
#SBATCH --nodes=1
##SBATCH --mem=200GB
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=save-order
#SBATCH --output=slurm-subfind-order.out
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
pigfile="/home1/08942/nianyic/ASTRID2_PIG/PIG_$snap"
dest="$dataroot/PIG_483_subfind"

ibrun -n 1 python3 ../code/1_save_subfind_order.py --pigfile $pigfile --subroot $subroot --tabfile $tabfile --grpfile $grpfile --dest $dest --cstart 338 --cend 341
