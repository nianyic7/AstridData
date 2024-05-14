#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 3:00:00
#SBATCH --job-name d2d-placeholder
#SBATCH --partition RM
#SBATCH -o slurm-%x.out
#SBATCH --mem=200G


source ~/.bashrc
conda activate default
codedir="$HOME/AstridData/bhdetails"

for snap in placeholder
do
echo "processing $snap"
ifile="/hildafs/datasets/Asterix/BH_details-chopped"
ofile="/hildafs/datasets/Asterix/BH_details_dict/first_pass"
srun python3 $codedir/detail_to_dict.py --idx $snap --srcdir $ifile --outdir $ofile 

ofile="/hildafs/datasets/Asterix/BH_details_dict/set2"
srun python3 $codedir/detail_to_dict2.py --idx $snap --srcdir $ifile --outdir $ofile

done
