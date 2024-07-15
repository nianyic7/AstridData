#!/bin/bash 
#SBATCH -p development
#SBATCH --time=2:00:00
#SBATCH --nodes=10
##SBATCH --mem=200GB
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nianyic@andrew.cmu.edu
#SBATCH --job-name=recover-021-gas
#SBATCH --output=%x.out

module load impi/19.0.9
export I_MPI_HYDRA_TOPOLIB=ipl
module unload python3
source activate fast-mpi4py

dataroot="$HOME/scratch3/Astrid/PART_temp"


snap="021"
ifile="$dataroot/PART_$snap-compressed"
ofile="$dataroot/PART_$snap"
codedir="$HOME/AstridData/compress_snapshot"

ibrun -n 280 python3 $codedir/fullsnap_recover.py --ifile $ifile --tol 0.03 --ofile $ofile --nfiles 256 --nloops 100 --blocknames "0/Position" "0/ElectronAbundance" "0/Mass" "0/InternalEnergy" "0/NeutralHydrogenFraction" "0/SmoothingLength"

wait
echo "Done with compressing"

