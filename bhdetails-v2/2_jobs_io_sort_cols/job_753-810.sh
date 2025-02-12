#!/bin/bash                                                                             
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=753,755,756,757,758,760,766,770,774,775,779,784,789,792,794,797,800,804,806,808,810
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
