#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time 2:00:00
#SBATCH --job-name chop-528-561
#SBATCH --partition development
#SBATCH -o slurm-%x.out
##SBATCH --mem=100G

source activate default
codedir="$HOME/AstridData/bhdetails"

for snap in 528 534 539 540 541 547 561
do
ifile="$HOME/scratch3/Astrid/bhdetails/BH-Details-R$snap"
ofile="$HOME/scratch3/Astrid/bhdetails-chopped/BH-Details-R$snap"
srun python3 $codedir/bhdetails_chop_full.py --ifile $ifile --ofile $ofile --nfiles 8 --nchunks 6
done

