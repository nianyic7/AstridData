#!/bin/bash                                                                             
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=493,496,497,500,506,508,509,514,518,523,528,534,539,540,541,547,552,557,561,566,572,577
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
