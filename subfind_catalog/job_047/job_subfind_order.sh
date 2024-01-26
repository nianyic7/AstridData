#!/bin/bash 
#SBATCH -p RM
#SBATCH --time=20:00:00
#SBATCH --nodes=2
#SBATCH --mem=100GB
#SBATCH --job-name=wseed
#SBATCH --output=slurm-subfind-order.out

snap=107
dataroot="/hildafs/datasets/Asterix"
tabfile="fof_subhalo_tab_$snap.hdf5"
grpfile="snap_$snap.hdf5"
subroot="$dataroot/subfind/subfind_$snap"
pigfile="$dataroot/PIG_files/PIG_$snap"
dest="$dataroot/PIG2/PIG_107_subfind"

mpirun -n 56 python3 ../save_subfind_order.py --pigfile $pigfile --subroot $subroot --tabfile $tabfile --grpfile $grpfile --dest $dest --cstart 0
