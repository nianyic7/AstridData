#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time 3:00:00
#SBATCH --job-name chop-placeholder
#SBATCH --partition RM
#SBATCH -o slurm-%x.out
#SBATCH --mem=200G

source ~/.bashrc
conda activate default

codedir="$HOME/AstridData/bhdetails"

for snap in placeholder
do
ifile="/hildafs/datasets/Asterix/BH_details_bigfile2/BH-Details-R$snap"
ofile="/hildafs/datasets/Asterix/BH_details-chopped/BH-Details-R$snap"
srun python3 $codedir/bhdetails_chop_full.py --ifile $ifile --ofile $ofile --nfiles 8 --nchunks 6
done

