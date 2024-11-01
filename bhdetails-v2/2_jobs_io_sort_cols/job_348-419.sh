#!/bin/bash                                                                             
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=348,353,359,364,369,370,374,375,376,378,381,386,388,392,395,399,405,408,414,419
#SBATCH --time 20:00:00
#SBATCH --job-name io-sort
#SBATCH --partition RM
#SBATCH -o slurm-%x-%a.out
#SBATCH --mem=200G

source ~/.bashrc
conda activate default

snap=$SLURM_ARRAY_TASK_ID
path="/hildafs/datasets/Asterix/BH_details_bigfile2"

srun  python3 ../io_sorted_cols.py --snap $snap --path $path
