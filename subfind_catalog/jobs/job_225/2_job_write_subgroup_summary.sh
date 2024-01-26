#!/bin/bash 
#SBATCH -p development
#SBATCH --time=2:00:00
#SBATCH --nodes=8
##SBATCH --mem=200GB
#SBATCH --ntasks-per-node=28
#SBATCH --job-name=summ
#SBATCH --output=slurm-write-summary.out

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

ibrun -n 224 python3 ../code/2_write_subgroup_summary.py --pigfile $pigfile --subroot $subroot --tabfile $tabfile --grpfile $grpfile --dest $dest --cstart 900 --minpart 200

#ibrun -n 224 python3 ../code/2_write_subgroup_summary.py --pigfile $pigfile --subroot $subroot --tabfile $tabfile --grpfile $grpfile --dest $dest --cstart 0 --minpart 20
