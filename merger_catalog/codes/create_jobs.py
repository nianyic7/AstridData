#!/usr/bin/python3
## file: submit_jobs.py
import os
import csv, subprocess
import numpy as np
import array
import struct
import argparse
from tempfile import mkstemp
from shutil import move, copymode
from os import fdopen, remove
import sys
import glob
import BigFile


################ write input files for AMBER ################
print("writing inputs...")
def replace(infile,outfile, pattern, subst):
    with open(outfile,'w') as new_file:
        with open(infile) as old_file:
            for line in old_file:
                for p, s in zip(pattern, subst):
                    line = line.replace(p, s)
                new_file.write(line)
        old_file.close()
    new_file.close()

    
def get_snaps_with_stars(roots):
    f4_snaps = []
    for root in roots:
        ffs = sorted(glob.glob(root + "PIG_*"))
        for ff in ffs:
            if os.path.exists(ff + "/4/Position") and os.path.exists(ff + "/5/Position"):
                pig       = BigFile(ff)
                battr     = pig["Header"].attrs
                scale_fac = pig["Header"].attrs["Time"][0]
                redshift  = 1./scale_fac - 1
                snap      = int(ff.split('_')[-1])

                f4_snaps.append((snap, redshift))
            
    dtype    = np.dtype([("snap_num", int), ("redshift", np.float32)])
    f4_snaps = np.array(f4_snaps, dtype=dtype)
    f4_snaps = np.sort(f4_snaps, order=["snap_num"])
    print("All snapshots with stars:", f4_snaps, flush=True)
    
    return f4_snaps


roots = ["/home1/08942/nianyic/ASTRID2_PIG/", "/home1/08942/nianyic/asterix/PIG/"]
allsnaps = get_snaps_with_stars(roots)

#---------------- Training Set -------------------------
# write input file 
samplefile = "inputs/job_merger_sample.sh"
root = "/home1/08942/nianyic/work/Astrid/merger_catalog/jobs/"



for ss in allsnaps:
    snap = ss[0]
    if snap <= 202:
        continue
    outfile = 'job_merger_%03d.sh'%(snap)
    full_out = root + outfile

     # write input files
    pattern = ['***snap***']
    subst = ['%d'%snap]
    for k in range(len(pattern)):
        print('Substitute ',pattern[k],'with',subst[k])
        
    replace(samplefile, full_out, pattern, subst)
    
    print('Finished writing:',outfile)
