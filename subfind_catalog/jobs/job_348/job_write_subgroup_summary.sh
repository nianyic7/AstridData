#!/bin/bash 
#SBATCH -p RM
#SBATCH --time=60:00:00
#SBATCH --nodes=4
#SBATCH --mem=200GB
#SBATCH --job-name=wseed3
#SBATCH --output=slurm-write-summary.out

snap=348
dataroot="/hildafs/datasets/Asterix"
tabfile="fof_subhalo_tab_$snap.hdf5"
grpfile="snap-groupordered_$snap.hdf5"
subroot="$dataroot/subfind/subfind_$snap"
pigfile="$dataroot/PIG2/PIG_$snap"
dest="$dataroot/PIG2/PIG_348_subfind"

mpirun -n 112 python3 ../write_subgroup_summary.py --pigfile $pigfile --subroot $subroot --tabfile $tabfile --grpfile $grpfile --dest $dest --cstart 800
