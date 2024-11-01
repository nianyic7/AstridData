#!/bin/bash                                                                             
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=5,19,29,38,45,46,52,55,58,64,68,72,73,75,83,86,88,92,97,99,102,103,110,114,119,125,127,131,135,140,141,146,150,154,157,161,162,164,169,175,180 
#SBATCH --time 20:00:00
#SBATCH --job-name io-sort
#SBATCH --partition RM
#SBATCH -o slurm-%x-%a.out
#SBATCH --mem=100G

source ~/.bashrc
conda activate default

snap=$SLURM_ARRAY_TASK_ID
path="/hildafs/datasets/Asterix/BH_details_bigfile"

srun  python3 ../io_sorted_cols.py --snap $snap --path $path
