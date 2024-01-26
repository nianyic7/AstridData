#!/bin/bash 
#SBATCH -p development
#SBATCH --time=2:00:00
#SBATCH --nodes 4
##SBATCH --mem=200GB
#SBATCH --ntasks-per-node=14
#SBATCH --job-name=wseed2
#SBATCH --output=slurm-group-reorder-4-5.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nianyic@andrew.cmu.edu

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

#ibrun -n 56 python3 ../code/3_group_reorder.py --pigfile $pigfile --dest $dest --gstart 0 --gend 30000 --blocknames "5/BlackholeAccretionRate" "5/BlackholeMass" "5/BlackholeMseed" "5/ID" "5/Position" "5/Velocity" "5/StarFormationTime"

#ibrun -n 56 python3 ../code/3_group_reorder.py --pigfile $pigfile --dest $dest --gstart 30000 --gend 300000 --blocknames "5/BlackholeAccretionRate" "5/BlackholeMass" "5/BlackholeMseed" "5/ID" "5/Position" "5/Velocity" "5/StarFormationTime"

#ibrun -n 56 python3 ../code/3_group_reorder.py --pigfile $pigfile --dest $dest --gstart 300000 --gend 3000000 --blocknames "5/BlackholeAccretionRate" "5/BlackholeMass" "5/BlackholeMseed" "5/ID" "5/Position" "5/Velocity" "5/StarFormationTime"

ibrun -n 56 python3 ../code/3_group_reorder.py --pigfile $pigfile --dest $dest --gstart 3000000 --gend 30000000 --blocknames "4/ID" "4/Mass" "4/Metallicity" "4/Metals" "4/Position" "4/StarFormationTime" "4/Velocity" "4/TotalMassReturned" "5/BlackholeAccretionRate" "5/BlackholeMass" "5/BlackholeMseed" "5/ID" "5/Position" "5/Velocity" "5/StarFormationTime"

