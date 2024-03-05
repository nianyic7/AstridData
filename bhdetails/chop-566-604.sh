#!/bin/bash
#SBATCH --nodes 1
#SBATCH -n 64
#SBATCH --time 4:00:00
#SBATCH --job-name chop-566-604
#SBATCH --partition RM-shared
#SBATCH -o slurm-%x.out
#SBATCH --mem=110G

source activate default
codedir="$HOME/AstridData/bhdetails"

for snap in 566 572 577 582 583 586 592 594 599 600 604
do
ifile="$HOME/scratch1/Astrid/bhdetails/BH-Details-R$snap"
ofile="$HOME/scratch1/Astrid/bhdetails-chopped/BH-Details-R$snap"
srun python3 $codedir/bhdetails_chop.py --ifile $ifile --ofile $ofile --nfiles 8 --nchunks 6
done

