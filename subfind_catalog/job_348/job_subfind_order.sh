#!/bin/bash 
#SBATCH -p RM
#SBATCH --time=20:00:00
#SBATCH --nodes=6
#SBATCH --mem=200GB
#SBATCH --job-name=wseed
#SBATCH --output=slurm-subfind-order.out

snap=348
dataroot="/hildafs/datasets/Asterix"
tabfile="fof_subhalo_tab_$snap.hdf5"
grpfile="snap-groupordered_$snap.hdf5"
subroot="$dataroot/subfind/subfind_$snap"
pigfile="$dataroot/PIG2/PIG_$snap"
dest="$dataroot/PIG2/PIG_348_subfind"

mpirun -n 168 python3 ../save_subfind_order.py --pigfile $pigfile --subroot $subroot --tabfile $tabfile --grpfile $grpfile --dest $dest --cstart 1
