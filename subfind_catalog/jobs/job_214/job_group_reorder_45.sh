#!/bin/bash 
#SBATCH -p RM
#SBATCH --time=48:00:00
#SBATCH --nodes 6
#SBATCH --mem=100GB
#SBATCH --job-name=wseed2
#SBATCH --output=slurm-group-reorder-4-5.out

snap=214
dataroot="/hildafs/datasets/Asterix"
tabfile="fof_subhalo_tab_$snap.hdf5"
grpfile="snap_$snap.hdf5"
subroot="$dataroot/subfind/subfind_$snap"
pigfile="$dataroot/PIG_files/PIG_$snap"
dest="$dataroot/PIG2/PIG_214_subfind"

mpirun -n 168 python3 ../code/3_group_reorder.py --gend 50000000 --pigfile $pigfile --dest $dest --blocknames "4/Mass" "4/Metallicity" "4/Position" "4/StarFormationTime" "4/Velocity" "5/BlackholeAccretionRate" "5/BlackholeMass" "5/BlackholeMseed" "5/ID" "5/Position" "5/Velocity" "5/StarFormationTime"
