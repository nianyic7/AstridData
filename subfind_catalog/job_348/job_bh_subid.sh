#!/bin/bash 
#SBATCH -p RM
#SBATCH --time=20:00:00
#SBATCH --nodes=4
#SBATCH --mem=100GB
#SBATCH --job-name=wseed3
#SBATCH --output=slurm-bh-sid.out

dataroot="/hildafs/datasets/Asterix"
dest="$dataroot/PIG2/PIG_348_subfind"

mpirun -n 112 python3 ../bh-subid.py  --dest $dest --gstart 0  --gend 30462029
