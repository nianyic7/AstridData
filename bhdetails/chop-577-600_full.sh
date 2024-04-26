#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=128
#SBATCH --time 2:00:00
#SBATCH --job-name chop-577-600-full
#SBATCH --partition development
#SBATCH -o slurm-%x.out
##SBATCH --mem=200G
source ~/.bashrc
conda activate default
codedir="$HOME/AstridData/bhdetails"

#for snap in 566 572 577 582 583 586 592 594 599 600 604
#do
#ifile="$HOME/scratch1/Astrid/bhdetails/BH-Details-R$snap"
#ofile="$HOME/scratch1/Astrid/bhdetails-chopped/BH-Details-R$snap"
#srun python3 $codedir/bhdetails_chop.py --ifile $ifile --ofile $ofile --nfiles 8 --nchunks 6
#done

for snap in 577 582 583 586 592 599 600
do
echo "processing $snap"
ifile="$HOME/scratch3/Astrid/bhdetail/BH-Details-R$snap"
ofile="$HOME/scratch3/Astrid/bhdetails-chopped/BH-Details-R$snap"
srun python3 $codedir/bhdetails_chop_full.py --ifile $ifile --ofile $ofile --nfiles 8 --nchunks 6
done
