#!/bin/bash 
#SBATCH -p RM
#SBATCH --time=20:00:00
#SBATCH --nodes=8
#SBATCH --mem=200GB
#SBATCH --job-name=wseed3
#SBATCH --output=slurm-append-offset.out

dataroot="/hildafs/datasets/Asterix"
dest="$dataroot/PIG2/PIG_348_subfind"

mpirun -n 224 python3 ../append_suboff.py  --dest $dest --gend 30462029
