#!/bin/bash 
#SBATCH -p development
#SBATCH --time=2:00:00
#SBATCH --nodes=8
##SBATCH --mem=200GB
#SBATCH --ntasks-per-node=28
#SBATCH --job-name=summ
#SBATCH --output=slurm-write-summary.out

module load impi/19.0.9
export I_MPI_HYDRA_TOPOLIB=ipl
module unload python3
source activate fast-mpi4py
which python

snap=505
tabfile="fof_subhalo_tab_$snap.hdf5"
grpfile="snap_$snap.hdf5"
subroot="/scratch3/09475/yihaoz/Astrid-subfind/$snap/subfind_$snap"
pigfile="/scratch3/09475/yihaoz/Astrid-subfind/$snap/PIG_$snap"
dest="$HOME/scratch3/Astrid/subfind-PIG/PIG_505_subfind"
codedir="$HOME/AstridData/subfind_catalog/code"


ibrun -n 224 python3 $codedir/2_write_subgroup_summary.py --pigfile $pigfile --subroot $subroot --tabfile $tabfile --grpfile $grpfile --dest $dest --cstart 1500 --minpart 200

#ibrun -n 224 python3 ../code/2_write_subgroup_summary.py --pigfile $pigfile --subroot $subroot --tabfile $tabfile --grpfile $grpfile --dest $dest --cstart 0 --minpart 20
