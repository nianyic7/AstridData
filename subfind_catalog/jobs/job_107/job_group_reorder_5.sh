#!/bin/bash 
#SBATCH -p RM
#SBATCH --time=48:00:00
#SBATCH --nodes 4   
#SBATCH --mem=100GB
#SBATCH --job-name=wseed2
#SBATCH --output=slurm-group-reorder-5.out

snap=107
dataroot="/hildafs/datasets/Asterix"
tabfile="fof_subhalo_tab_$snap.hdf5"
grpfile="snap_$snap.hdf5"
subroot="$dataroot/subfind/subfind_$snap"
pigfile="$dataroot/PIG_files/PIG_$snap"
dest="$dataroot/PIG2/PIG_107_subfind"



mpirun -n 112 python3 ../group_reorder.py --pigfile $pigfile --dest $dest --minpart 1 --gstart 0 --gend 1000000000 --blocknames  "5/BlackholeAccretionRate" "5/BlackholeMass" "5/BlackholeMseed" "5/ID" "5/Position" "5/Velocity" "5/StarFormationTime"
