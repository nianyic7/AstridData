#!/bin/bash 
#SBATCH -p RM
#SBATCH --time=48:00:00
#SBATCH --nodes 12
#SBATCH --mem=200GB
#SBATCH --job-name=wseed2
#SBATCH --output=slurm-group-reorder-1.out

snap=214
dataroot="/hildafs/datasets/Asterix"
tabfile="fof_subhalo_tab_$snap.hdf5"
grpfile="snap_$snap.hdf5"
subroot="$dataroot/subfind/subfind_$snap"
pigfile="$dataroot/PIG2/PIG_$snap"
dest="$dataroot/PIG2/PIG_214_subfind"


mpirun -n 336 python3 ../code/3_group_reorder.py --pigfile $pigfile --dest $dest --gstart 0 --gend 50000 --blocknames  "1/Position" "1/Velocity"
mpirun -n 336 python3 ../code/3_group_reorder.py --pigfile $pigfile --dest $dest --gstart 50000 --gend 200000 --blocknames  "1/Position" "1/Velocity"
mpirun -n 336 python3 ../code/3_group_reorder.py --pigfile $pigfile --dest $dest --gstart 200000 --gend 2000000 --blocknames "1/Position" "1/Velocity"
mpirun -n 336 python3 ../code/3_group_reorder.py --pigfile $pigfile --dest $dest --gstart 2000000 --gend 8000000 --blocknames  "1/Position" "1/Velocity"
mpirun -n 336 python3 ../code/3_group_reorder.py --pigfile $pigfile --dest $dest --gstart 8000000 --gend 50000000 --blocknames  "1/Position" "1/Velocity"
