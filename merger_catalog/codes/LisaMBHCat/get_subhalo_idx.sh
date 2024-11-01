#!/bin/bash                                                                                                                                                                                                                                                     
#SBATCH -p RM
#SBATCH -N 2
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=7
#SBATCH --job-name=sidx
#SBATCH --time=10:00:00
#SBATCH --mem=200GB
#SBATCH --output=slurm-save-sidx.log

#source ~/.bashrc
#conda activate default
subroot="/hildafs/datasets/Asterix/PIG2"
mfile="/hildafs/home/nianyic/scratch1/Astrid_data/lisa_mbhcat/up_to_date_mergers_nogal.npy"
outdir="/hildafs/home/nianyic/scratch1/Astrid_data/lisa_mbhcat/sidx"

mpirun -n 8 python merger_hostidx.py --subroot $subroot  --mfile $mfile  --outdir $outdir

