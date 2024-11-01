#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 8:00:00
#SBATCH --job-name d2d-514
#SBATCH --partition HENON
#SBATCH -o slurm-%x.out
#SBATCH --mem=400G


source ~/.bashrc
conda activate default
codedir="$HOME/AstridData/bhdetails"

for snap in 514
do
echo "processing $snap"
ifile="/hildafs/datasets/Asterix/BH_details_bigfile2"
ofile="/hildafs/datasets/Asterix/BH_details_484-653/first_pass_minimal"
srun python3 $codedir/detail_to_dict_minimal.py --idx $snap --srcdir $ifile --outdir $ofile 

done
