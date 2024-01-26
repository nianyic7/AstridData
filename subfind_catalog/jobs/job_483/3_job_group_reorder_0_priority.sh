#!/bin/bash 
#SBATCH -p development
#SBATCH --time=2:00:00
#SBATCH --nodes 16
##SBATCH --mem=200GB
#SBATCH --ntasks-per-node=14
#SBATCH --job-name=wseed2
#SBATCH --output=slurm-group-reorder-0.out
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

#ibrun -n 224 python3 ../code/3_group_reorder.py --pigfile $pigfile --dest $dest --gstart 0 --gend 30000 --blocknames "0/Position"  "0/Mass" "0/StarFormationRate" "0/SmoothingLength" "0/Metallicity" "0/InternalEnergy" "0/NeutralHydrogenFraction" "0/Density" "0/DelayTime" "0/ElectronAbundance"
sleep 1
echo "Finished First Batch"

#ibrun -n 224 python3 ../code/3_group_reorder.py --pigfile $pigfile --dest $dest --gstart 30000 --gend 300000 --blocknames "0/Position"  "0/Mass" "0/StarFormationRate" "0/SmoothingLength" "0/Metallicity" "0/InternalEnergy" "0/NeutralHydrogenFraction" "0/Density" "0/DelayTime" "0/ElectronAbundance"
sleep 1
echo "Finished Second Batch"

ibrun -n 224 python3 ../code/3_group_reorder.py --pigfile $pigfile --dest $dest --gstart 300000 --gend 1000000 --blocknames "0/Position"  "0/Mass" "0/StarFormationRate" "0/SmoothingLength" "0/Metallicity" "0/InternalEnergy" "0/NeutralHydrogenFraction" "0/Density" "0/DelayTime" "0/ElectronAbundance"
sleep 1
echo "Finished Third Batch"

ibrun -n 224 python3 ../code/3_group_reorder.py --pigfile $pigfile --dest $dest --gstart 1000000 --gend 3000000 --blocknames "0/Position"  "0/Mass" "0/StarFormationRate" "0/SmoothingLength" "0/Metallicity" "0/InternalEnergy" "0/NeutralHydrogenFraction" "0/Density" "0/DelayTime" "0/ElectronAbundance"
sleep 1
echo "Finished Fourth Batch"


#ibrun -n 224 python3 ../code/3_group_reorder.py --pigfile $pigfile --dest $dest --gstart 3000000 --gend 30000000 --blocknames "0/Position"  "0/Mass" "0/StarFormationRate" "0/SmoothingLength" "0/Metallicity" "0/InternalEnergy" "0/NeutralHydrogenFraction" "0/Density" "0/DelayTime" "0/ElectronAbundance"
sleep 1
echo "Finished Last Batch"
