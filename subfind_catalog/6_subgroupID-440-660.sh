#!/bin/bash 
#SBATCH -p RM
#SBATCH --time=2:00:00
#SBATCH --nodes=1
##SBATCH --mem=200GB
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=subidx-440-660
#SBATCH --output=slurm-%x.out
#SBATCH --mail-type=END
#SBATCH --mail-user=nianyic@andrew.cmu.edu


module load intelmpi gsl intel
export I_MPI_HYDRA_TOPOLIB=ipl

source ~/.bashrc
conda activate default
which python


particle_type=5
dataroot="/hildafs/datasets/Asterix/PIG2/"
codedir="/hildafs/home/nianyic/AstridData/subfind_catalog/code"

for snap in 440 461 483 505 523 544 567 590 609 634 658 660 
do

    indir="$dataroot/PIG_${snap}_subfind"
    odir="$dataroot/SubIdx/S-${snap}"

    # if the problematic SubgroupIndex was generated used, mv it to SubgroupIndex_old

    if [ -d $indir/$particle_type/SubgroupIndex_test ]; then
        mv $indir/$particle_type/SubgroupIndex_test $indir/$particle_type/SubgroupIndex_old
        echo "moving $indir/$particle_type/SubgroupIndex_test to $indir/$particle_type/SubgroupIndex_old"
    fi

    python3 $codedir/6_bh-subid.py  --src $indir --dest $odir --particle_type $particle_type 

done

