#!/bin/bash
source activate default
codedir="$HOME/AstridData/bhdetails"

for snap in 528 534 539 540 541 547 561
do
echo "processing $snap"
ifile="$HOME/scratch1/Astrid/bhdetails/BH-Details-R$snap"
ofile="$HOME/scratch1/Astrid/bhdetails-chopped/BH-Details-R$snap"
python3 $codedir/bhdetails_chop.py --ifile $ifile --ofile $ofile --nfiles 32 --nchunks 6
done

