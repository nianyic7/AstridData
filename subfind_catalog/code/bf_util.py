import numpy as np
from bigfile import BigFile
import h5py
import sys,os
import glob
import argparse



def dtype_size_nfile(pig,blockname):
    return pig[blockname].dtype, pig[blockname].size, pig[blockname].Nfile


def init_block(dest,blockname,dtype,size,nfile):
    try:
        block = dest[blockname]
    except:
        dest.create(blockname, dtype=dtype, size=size, Nfile=nfile)
        block = dest[blockname]
    return block
        
def get_features(pigdir,ptype):
    features = []
    for f in glob.glob(pigdir + '/%d/*'%ptype):
        features.append(f.split('/')[-1])
    return features
    
    





