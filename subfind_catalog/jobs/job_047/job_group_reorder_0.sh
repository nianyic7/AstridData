#!/bin/bash 
#SBATCH -p RM
#SBATCH --time=48:00:00
#SBATCH --nodes 6   
#SBATCH --mem=200GB
#SBATCH --job-name=wseed2
#SBATCH --output=slurm-group-reorder-0.out

snap=047
dataroot="/hildafs/datasets/Asterix"
tabfile="fof_subhalo_tab_$snap.hdf5"
grpfile="snap_$snap.hdf5"
subroot="$dataroot/subfind/subfind_$snap"
pigfile="$dataroot/PIG_files/PIG_$snap"
dest="$dataroot/PIG2/PIG_047_subfind"


mpirun -n 168 python3 ../group_reorder.py --pigfile $pigfile --dest $dest --minpart 1 --gstart 0 --gend 50000 --blocknames "0/Position" "0/Velocity" "0/Mass" "0/StarFormationRate" "0/SmoothingLength" "0/Metallicity" "0/InternalEnergy" "0/NeutralHydrogenFraction" "0/Density" "0/EgyWtDensity" "0/ElectronAbundance" "0/ID"
mpirun -n 168 python3 ../group_reorder.py --pigfile $pigfile --dest $dest --minpart 1 --gstart 50000 --gend 200000 --blocknames  "0/Position" "0/Velocity" "0/Mass" "0/StarFormationRate" "0/SmoothingLength" "0/Metallicity" "0/InternalEnergy" "0/NeutralHydrogenFraction" "0/Density" "0/EgyWtDensity" "0/ElectronAbundance" "0/ID" 

mpirun -n 168 python3 ../group_reorder.py --pigfile $pigfile --dest $dest --minpart 1 --gstart 200000 --gend 2000000 --blocknames  "0/Position" "0/Velocity" "0/Mass" "0/StarFormationRate" "0/SmoothingLength" "0/Metallicity" "0/InternalEnergy" "0/NeutralHydrogenFraction" "0/Density" "0/EgyWtDensity" "0/ElectronAbundance" "0/ID" 

mpirun -n 168 python3 ../group_reorder.py --pigfile $pigfile --dest $dest --minpart 1 --gstart 2000000 --gend 8000000 --blocknames  "0/Position" "0/Velocity" "0/Mass" "0/StarFormationRate" "0/SmoothingLength" "0/Metallicity" "0/InternalEnergy" "0/NeutralHydrogenFraction" "0/Density" "0/EgyWtDensity" "0/ElectronAbundance" "0/ID" 

mpirun -n 168 python3 ../group_reorder.py --pigfile $pigfile --dest $dest --minpart 1 --gstart 8000000 --gend 50000000 --blocknames  "0/Position" "0/Velocity" "0/Mass" "0/StarFormationRate" "0/SmoothingLength" "0/Metallicity" "0/InternalEnergy" "0/NeutralHydrogenFraction" "0/Density" "0/EgyWtDensity" "0/ElectronAbundance" "0/ID" 

mpirun -n 168 python3 ../group_reorder.py --pigfile $pigfile --dest $dest --minpart 1 --gstart 50000000 --gend 1867081150 --blocknames  "0/Position" "0/Velocity" "0/Mass" "0/StarFormationRate" "0/SmoothingLength" "0/Metallicity" "0/InternalEnergy" "0/NeutralHydrogenFraction" "0/Density" "0/EgyWtDensity" "0/ElectronAbundance" "0/ID" 




