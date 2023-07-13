#!/bin/bash 
#SBATCH -p development
#SBATCH --time=2:00:00
#SBATCH --nodes=20
##SBATCH --mem=200GB
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nianyic@andrew.cmu.edu
#SBATCH --job-name=comp494
#SBATCH --output=slurm-%x.out

module load impi/19.0.9
export I_MPI_HYDRA_TOPOLIB=ipl
module unload python3
source activate fast-mpi4py

dataroot="$HOME/ASTRID2_PIG"
outroot="$HOME/scratch3/Astrid/PART_compressed"
mkdir -p $outroot


for snap in 494 
do
ifile="$dataroot/PART_$snap"
ofile="$outroot/PART_$snap-compressed"
codedir="$HOME/AstridData/compression/"

ibrun -n 560 python3 $codedir/fullsnap_compress.py --ifile $ifile --tol 0.03 --ofile $ofile --nfiles 256 --nloops 100 --blocknames  "1/Position" "1/Velocity" "0/Position" "4/Position"
done

wait
echo "Done with compressing"

date
