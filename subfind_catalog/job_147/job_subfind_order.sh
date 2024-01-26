#!/bin/bash 
#SBATCH -p RM
#SBATCH --time=20:00:00
#SBATCH --nodes=4
#SBATCH --mem=100GB
#SBATCH --job-name=sf
#SBATCH --output=slurm-subfind-order.out


snap=147
dataroot="/hildafs/datasets/Asterix"
tabfile="fof_subhalo_tab_$snap.hdf5"
grpfile="snap_$snap.hdf5"
subroot="$dataroot/subfind/subfind_$snap"
pigfile="$dataroot/PIG_files/PIG_$snap"
dest="$dataroot/PIG2/PIG_147_subfind"

mpirun -n 112 python3 ../script147/save_subfind_order.py --snap $snap --pigfile $pigfile --subroot $subroot --tabfile $tabfile --grpfile $grpfile --dest $dest --cstart 0
