import numpy as np
from bigfile import BigFile
import struct
import glob
import os, sys
from mpi4py import MPI
from bigfile import FileMPI
import argparse
from natsort import natsorted



SEL_COL = ['BHID','BHMass','Mdot','Density','timebin','Encounter','MinPos',\
             'MinPot','Entropy','GasVel','acMom','acMass','acBHMass',\
             'Fdbk','SPHID','SwallowID','CountProgs','Swallowed',\
             'BHpos','srDensity','srParticles','srVel','srDisp',\
             'DFAccel','DragAccel','GravAccel','BHvel','Mtrack','Mdyn',\
             'KineticFdbkEnergy','NumDM','V1sumDM','V2sumDM','MgasEnc','KEflag','z']


class BHType:
    def __init__(self, name_sel=None):
        self.name_all = ('size','BHID','BHMass','Mdot','Density','timebin','Encounter','MinPos',\
             'MinPot','Entropy','GasVel','acMom','acMass','acBHMass',\
             'Fdbk','SPHID','SwallowID','CountProgs','Swallowed',\
             'BHpos','srDensity','srParticles','srVel','srDisp',\
             'DFAccel','DragAccel','GravAccel','BHvel','Mtrack','Mdyn',\
             'KineticFdbkEnergy','NumDM','V1sumDM','V2sumDM','MgasEnc','KEflag','z','size2')
        self.dtype_all = ('i','q','d','d','d','i','i','3d',\
            'd','d','3d','3d','d','d',\
            'd','q','q','i','i',\
            '3d','d','d','3d','d',\
            '3d','3d','3d','3d','d','d',\
            'd','d','3d','d','d','i','d','i')
        self.dtype_lean = ('i','q','f','f','f','i','i','3f',\
            'f','f','3f','3f','f','f',\
            'f','q','q','i','i',\
            '3f','f','f','3f','f',\
            '3f','3f','3f','3f','f','f',\
            'f','f','3f','f','f','i','f','i')

        if name_sel is None:
            name_sel = self.name_all
        self.name_sel = name_sel
        self.name2type = {name: dtype for name, dtype in zip(self.name_all, self.dtype_all)}
        self.name2type_lean = {name: dtype for name, dtype in zip(self.name_all, self.dtype_lean)}


    @property
    def TypeAll(self):
        np_type = np.dtype({'names':self.name_all, 'formats':self.dtype_all})
        return np_type
    
    @property
    def TypeLean(self):
        np_type = np.dtype({'names':self.name_all, 'formats':self.dtype_lean})
        return np_type
    
    @property
    def TypeSel(self):
        name_sel = [name for name in self.name_all if name in name_sel]
        dtype_sel = [dtype for name, dtype in zip(self.name_all, self.dtype_all) if name in name_sel]
        np_type = np.dtype({'names':name_sel, 'formats':dtype_sel})
        return np_type


def init_block(bf_w, blocknames, nfile):
    block_dict = {}
    for blockname in blocknames:
        dtype = BHTYPE.name2type_lean[blockname]
        dsize = NTOT
        block = bf_w.create(blockname, dtype, dsize, nfile)
        block_dict[blockname] = block
        print("Initizlized block:", blockname, flush=True)
    return block_dict

if __name__ == "__main__":

    # -------------- Cmd line Args ------------------------------
    
    parser = argparse.ArgumentParser(description="snap")
    parser.add_argument(
        "--fileroot",
        required=True,
        type=str,
        help="path of the input file (end without /)",
    )
    parser.add_argument(
        "--snap", required=True, type=int, help="snapshot number"
    )
    parser.add_argument(
        "--ofile",
        required=True,
        type=str,
        help="path of the output file directory (end without /)",
    )
    parser.add_argument(
        "--nfiles", required=True, type=int, help="number of files in the output snap"
    )

    args = parser.parse_args()

    # --------------------------
    BHTYPE = BHType(SEL_COL)

    ifilelist = natsorted(glob.glob(f"{args.fileroot}/BlackholeDetails-R{args.snap:03d}*"))
    if len(ifilelist) == 0:
        print(f"No file to process", flush=True)
        sys.exit(0)

    # create ofile
    bf_w = BigFile(args.ofile, create=True)
    # get total data size
    len_list = []
    off_list = []
    offset = 0
    for file in ifilelist:
        print(f"Reading file: {file}", flush=True)
        bf = BigFile(file)
        size = bf["BHID"].size
        print(f"Size of file: {size}", flush=True)
        len_list.append(size)
        off_list.append(offset)
        offset += len_list[-1]
    NTOT = sum(len_list)


    # initialize block
    block_dict = init_block(bf_w, blocknames=SEL_COL, nfile=args.nfiles)
    # io block
    for i, file in enumerate(ifilelist):
        bf = BigFile(file)
        print(f"Processing file {file}", flush=True)
        for blockname in SEL_COL:
            print(f"Block {blockname}", flush=True)
            block = block_dict[blockname]
            block.write(off_list[i], bf[blockname][:])
        print(f"Finished processing file {file}", flush=True)
    
    print("Done!", flush=True)
