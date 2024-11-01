#!/bin/bash                                                                                                                                                                                                                                                     
#SBATCH -p RM
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --job-name=sidx
#SBATCH --time=10:00:00
#SBATCH --mem=100GB
#SBATCH --output=slurm-test-214.log 
source ~/.bashrc
conda activate default
subroot="/hildafs/datasets/Asterix/PIG2"
mfile="/hildafs/home/nianyic/scratch1/Astrid_data/lisa_mbhcat/up_to_date_mergers_nogal.npy"
outdir="/hildafs/home/nianyic/scratch1/Astrid_data/lisa_mbhcat/sidx"
 
srun python3 merger_hostidx_single_snap.py --snap 214 --subroot $subroot  --mfile $mfile  --outdir $outdir
