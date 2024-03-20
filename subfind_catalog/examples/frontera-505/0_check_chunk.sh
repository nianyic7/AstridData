#!/bin/bash 
#SBATCH -p development
#SBATCH --time=1:00:00
#SBATCH --nodes=2
##SBATCH --mem=200GB
#SBATCH --ntasks-per-node=14
#SBATCH --job-name=check-chunk
#SBATCH --output=slurm-check-chunk.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nianyic@andrew.cmu.edu

module unload python3
source activate fast-mpi4py 
which python
snap=225
dataroot="/home1/08942/nianyic/scratch3/Astrid"
tabfile="fof_subhalo_tab_$snap.hdf5"
grpfile="snap_$snap.hdf5"
subroot="$dataroot/subfind/subfind-$snap"
savedir="$dataroot/subfind/checks-$snap"
codedir="$HOME/AstridData/subfind_catalog/code"
cd $savedir
rm *
ibrun -n 28 python3 $codedir/0_check_chunk.py --subroot $subroot --tabfile $tabfile --grpfile $grpfile --savedir $savedir --minsmass 2e6 
