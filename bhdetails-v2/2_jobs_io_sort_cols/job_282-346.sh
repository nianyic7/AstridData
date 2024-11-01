#!/bin/bash                                                                             
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=282,284,288,292,297,302,303,304,308,312,316,318,323,327,332,333,335,341,346
#SBATCH --time 20:00:00
#SBATCH --job-name io-sort
#SBATCH --partition RM
#SBATCH -o slurm-%x-%a.out
#SBATCH --mem=200G

source ~/.bashrc
conda activate default

snap=$SLURM_ARRAY_TASK_ID
path="/hildafs/datasets/Asterix/BH_details_bigfile"

srun  python3 ../io_sorted_cols.py --snap $snap --path $path
