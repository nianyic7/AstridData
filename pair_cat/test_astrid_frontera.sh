#!/bin/bash
module unload python3
source activate fast-mpi4py

python3 save_astrid_cat_hdf5.py --cluster "frontera" --snap 557 --Lmin 0 --Mmin 1e7 --rmax 5 --odir /home1/08942/nianyic/work/Astrid/MBHPair_catalog 

