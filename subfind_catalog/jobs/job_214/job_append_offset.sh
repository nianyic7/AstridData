#!/bin/bash 
#SBATCH -p RM
#SBATCH --time=20:00:00
#SBATCH --nodes=4
#SBATCH --mem=200GB
#SBATCH --job-name=wseed3
#SBATCH --output=slurm-append-offset.out

dataroot="/hildafs/datasets/Asterix"
dest="$dataroot/PIG2/PIG_214_subfind"

mpirun -n 112 python3 ../append_suboff.py  --dest $dest --gend 35999999
