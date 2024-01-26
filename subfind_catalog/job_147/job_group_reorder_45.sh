#!/bin/bash 
#SBATCH -p RM
#SBATCH --time=48:00:00
#SBATCH --nodes 8
#SBATCH --mem=200GB
#SBATCH --job-name=wseed2
#SBATCH --output=slurm-group-reorder-45.out

snap=147
dataroot="/hildafs/datasets/Asterix"
tabfile="fof_subhalo_tab_$snap.hdf5"
grpfile="snap_$snap.hdf5"
subroot="$dataroot/subfind/subfind_$snap"
pigfile="$dataroot/PIG_files/PIG_$snap"
dest="$dataroot/PIG2/PIG_147_subfind"



mpirun -n 224 python3 ../group_reorder.py --pigfile $pigfile --dest $dest --minpart 1 --gstart 0 --gend 1000000000 --blocknames "4/Mass" "4/Metallicity" "4/Position" "4/StarFormationTime" "4/Velocity" "4/ID" "5/BlackholeAccretionRate" "5/BlackholeMass" "5/BlackholeMseed" "5/ID" "5/Position" "5/Velocity" "5/StarFormationTime"
