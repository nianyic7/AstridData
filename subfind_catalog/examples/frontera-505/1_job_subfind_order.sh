#!/bin/bash 
#SBATCH -p development
#SBATCH --time=2:00:00
#SBATCH --nodes=4
##SBATCH --mem=200GB
#SBATCH --ntasks-per-node=14
#SBATCH --job-name=save-order
#SBATCH --output=slurm-subfind-order.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nianyic@andrew.cmu.edu

module load impi/19.0.9
export I_MPI_HYDRA_TOPOLIB=ipl
module unload python3
source activate fast-mpi4py

snap=505
tabfile="fof_subhalo_tab_$snap.hdf5"
grpfile="snap_$snap.hdf5"
subroot="/scratch3/09475/yihaoz/Astrid-subfind/$snap/subfind_$snap"
pigfile="/scratch3/09475/yihaoz/Astrid-subfind/$snap/PIG_$snap"
dest="$HOME/scratch3/Astrid/subfind-PIG/PIG_505_subfind"

codedir="$HOME/AstridData/subfind_catalog/code"



ibrun -n 56 python3 $codedir/1_save_subfind_order.py --pigfile $pigfile --subroot $subroot --tabfile $tabfile --grpfile $grpfile --dest $dest --cstart 1000 --cend 1200
