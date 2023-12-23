import os
import glob
import csv, subprocess
import numpy as np
import array
import struct
import argparse
from tempfile import mkstemp
from shutil import move, copymode
from os import fdopen, remove


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
    
    
root =  '/hildafs/datasets/Asterix/PIG_files/'
for subdir in sorted(glob.glob(root+'*/')):
    snap = subdir[-4:-1]
    if int(snap)<=147:
        continue
    
    samplejob = 'job_template.sh'
    jobfile = 'job-'+snap + '.sh'
    pattern = ['sub-snap']
    subst = [snap]
    replace(samplejob,jobfile, pattern, subst)
    print('Finished writing:',jobfile)

        
        
            
        