#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time 2:00:00
#SBATCH --job-name chop-642-653
#SBATCH --partition development
#SBATCH -o slurm-%x.out
##SBATCH --mem=200G

source activate default
codedir="$HOME/AstridData/bhdetails"

for snap in 642 647 648 653
do
ifile="$HOME/scratch3/Astrid/bhdetails/BH-Details-R$snap"
ofile="$HOME/scratch3/Astrid/bhdetails-chopped/BH-Details-R$snap"
srun python3 $codedir/bhdetails_chop_full.py --ifile $ifile --ofile $ofile --nfiles 8 --nchunks 6
done

