#!/bin/bash 
#SBATCH -p development
#SBATCH --time=2:00:00
#SBATCH --nodes=20
##SBATCH --mem=200GB
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nianyic@andrew.cmu.edu
#SBATCH --job-name=comp-489
#SBATCH --output=%x.out

module load impi/19.0.9
export I_MPI_HYDRA_TOPOLIB=ipl
module unload python3
source activate fast-mpi4py

outroot="$HOME/scratch3/Astrid/PIG_compressed"
dataroot="$HOME/ASTRID2_PIG"


snap=489
ifile="$dataroot/PIG_$snap"
ofile="$outroot/PIG_$snap-compressed"
codedir="$HOME/AstridData/compress_snapshot"

ibrun -n 560 python3 $codedir/fullsnap_compress.py --ifile $ifile --tol 0.03 --ofile $ofile --nfiles 256 --nloops 100

wait
echo "Done with compressing"

cp -r $ifile/Header $ofile
cp -r $ifile/FOFGroups $ofile

echo "copied over header and FOF"
date
