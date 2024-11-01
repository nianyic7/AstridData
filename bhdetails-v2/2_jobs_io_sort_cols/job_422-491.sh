#!/bin/bash                                                                             
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=422,427,433,434,438,444,447,453,458,464,469,475,477,478,479,481,484,486,491
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
