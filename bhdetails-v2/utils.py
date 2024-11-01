import numpy as np
from bigfile import BigFile
import struct
import glob
import os, sys

def get_filenum_list():
    path1 = "/hildafs/datasets/Asterix/BH_details_bigfile"
    path2 = "/hildafs/datasets/Asterix/BH_details_bigfile2"
    flist1 = sorted(glob.glob(path1 + "/BH-Details-R*"))
    flist2 = sorted(glob.glob(path2 + "/BH-Details-R*"))
    slist1 = [int(f.split("-")[-1][-3:]) for f in flist1]
    slist2 = [int(f.split("-")[-1][-3:]) for f in flist2]
    return slist1 + slist2

print(get_filenum_list())
