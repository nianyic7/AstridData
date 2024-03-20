#!/bin/bash 
#SBATCH -p development
#SBATCH --time=2:00:00
#SBATCH --nodes 16
##SBATCH --mem=200GB
#SBATCH --ntasks-per-node=14
#SBATCH --job-name=group-order-01
#SBATCH --output=slurm-group-reorder-0-1.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nianyic@andrew.cmu.edu

module load impi/19.0.9
export I_MPI_HYDRA_TOPOLIB=ipl
module unload python3
source activate fast-mpi4py
which python

snap=505
tabfile="fof_subhalo_tab_$snap.hdf5"
grpfile="snap_$snap.hdf5"
subroot="/scratch3/09475/yihaoz/Astrid-subfind/$snap/subfind_$snap"
pigfile="/scratch3/09475/yihaoz/Astrid-subfind/$snap/PIG_$snap"
dest="$HOME/scratch3/Astrid/subfind-PIG/PIG_505_subfind"
codedir="$HOME/AstridData/subfind_catalog/code"


# ignoring Generation
ibrun -n 224 python3 $codedir/3_group_reorder.py --pigfile $pigfile --dest $dest --gstart 0 --gend 30000 --blocknames "0/Position"  "0/Mass" "0/StarFormationRate" "0/SmoothingLength" "0/Metallicity" "0/InternalEnergy" "0/NeutralHydrogenFraction" "0/Density" "0/DelayTime" "0/ElectronAbundance" "0/EgyWtDensity" "0/Velocity" "0/HeIIIIonized" "0/Metals" "0/ID" "0/GroupID"
sleep 1
echo "Finished First Batch"

ibrun -n 224 python3 $codedir/3_group_reorder.py --pigfile $pigfile --dest $dest --gstart 30000 --gend 300000 --blocknames "0/Position"  "0/Mass" "0/StarFormationRate" "0/SmoothingLength" "0/Metallicity" "0/InternalEnergy" "0/NeutralHydrogenFraction" "0/Density" "0/DelayTime" "0/ElectronAbundance" "0/EgyWtDensity" "0/Velocity" "0/HeIIIIonized" "0/Metals" "0/ID" "0/GroupID"
sleep 1
echo "Finished Second Batch"

ibrun -n 224 python3 $codedir/3_group_reorder.py --pigfile $pigfile --dest $dest --gstart 300000 --gend 1000000 --blocknames "0/Position"  "0/Mass" "0/StarFormationRate" "0/SmoothingLength" "0/Metallicity" "0/InternalEnergy" "0/NeutralHydrogenFraction" "0/Density" "0/DelayTime" "0/ElectronAbundance" "0/EgyWtDensity" "0/Velocity" "0/HeIIIIonized" "0/Metals" "0/ID" "0/GroupID"
sleep 1
echo "Finished Third Batch"

#ibrun -n 224 python3 ../code/3_group_reorder.py --pigfile $pigfile --dest $dest --gstart 1000000 --gend 3000000 --blocknames "0/Position"  "0/Mass" "0/StarFormationRate" "0/SmoothingLength" "0/Metallicity" "0/InternalEnergy" "0/NeutralHydrogenFraction" "0/Density" "0/DelayTime" "0/ElectronAbundance" "0/EgyWtDensity" "0/Velocity" "0/HeIIIIonized" "0/Metals" "0/ID" "0/GroupID"
#sleep 1
#echo "Finished Fourth Batch"


#ibrun -n 224 python3 ../code/3_group_reorder.py --pigfile $pigfile --dest $dest --gstart 3000000 --gend 30000000 --blocknames "0/Position"  "0/Mass" "0/StarFormationRate" "0/SmoothingLength" "0/Metallicity" "0/InternalEnergy" "0/NeutralHydrogenFraction" "0/Density" "0/DelayTime" "0/ElectronAbundance" "0/EgyWtDensity" "0/Velocity" "0/HeIIIIonized" "0/Metals" "0/ID" "0/GroupID"
#sleep 1
#echo "Finished Last Batch"
