#!/bin/bash 
#SBATCH -p RM
#SBATCH --time=20:00:00
#SBATCH --nodes=6
#SBATCH --mem=200GB
#SBATCH --job-name=wseed3
#SBATCH --output=slurm-write-summary.out

snap=107
dataroot="/hildafs/datasets/Asterix"
tabfile="fof_subhalo_tab_$snap.hdf5"
grpfile="snap_$snap.hdf5"
subroot="$dataroot/subfind/subfind_$snap"
pigfile="$dataroot/PIG_files/PIG_$snap"
dest="$dataroot/PIG2/PIG_107_subfind"

mpirun -n 168 python3 ../write_subgroup_summary.py --pigfile $pigfile --subroot $subroot --tabfile $tabfile --grpfile $grpfile --dest $dest
