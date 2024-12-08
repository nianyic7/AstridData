#!/bin/bash                                                                             
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=734,737,738,741,744,747,749,751 
#SBATCH --time 20:00:00
#SBATCH --job-name io-sort
#SBATCH --partition HENON
#SBATCH -o slurm-%x-%a.out
#SBATCH --mem=200G

source ~/.bashrc
conda activate default

snap=$SLURM_ARRAY_TASK_ID
path="/hildafs/datasets/Asterix/BH_details_bigfile2"

srun  python3 ../2_io_sorted_cols.py --snap $snap --path $path
