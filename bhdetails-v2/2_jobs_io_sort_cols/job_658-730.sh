#!/bin/bash                                                                             
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=658,660,665,670,674,678,683,688,694,695,697,698,703,706,710,714,716,720,723,725,728,730 
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
