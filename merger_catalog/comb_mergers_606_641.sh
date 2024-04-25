#!/bin/bash
idir="/jet/home/nianyic/scratch1/Astrid/bhdetails-chopped/mergers"
odir="/jet/home/nianyic/scratch1/Astrid/mergers"

python combine_chopped_mergers.py --mdir $idir --smin 606 --smax 641 --outdir $odir --plot 1
