import numpy as np
from bigfile import BigFile
import struct
import glob
import os, sys
import argparse

name_all = ('size','BHID','BHMass','Mdot','Density','timebin','Encounter','MinPos',\
        'MinPot','Entropy','GasVel','acMom','acMass','acBHMass',\
        'Fdbk','SPHID','SwallowID','CountProgs','Swallowed',\
        'BHpos','srDensity','srParticles','srVel','srDisp',\
        'DFAccel','DragAccel','GravAccel','BHvel','Mtrack','Mdyn',\
        'KineticFdbkEnergy','NumDM','V1sumDM','V2sumDM','MgasEnc','KEflag','z','size2')
dtype_all = ('i','q','d','d','d','i','i','3d',\
    'd','d','3d','3d','d','d',\
    'd','q','q','i','i',\
    '3d','d','d','3d','d',\
    '3d','3d','3d','3d','d','d',\
    'd','d','3d','d','d','i','d','i')
name2type = {name: dtype for name, dtype in zip(name_all, dtype_all)}

SEL_COL=list(name_all)

def get_sort_indices(redshift):
    return np.argsort(redshift)[::-1]


def get_split_indices(redshift, sort_indices, nchunks):
    n = len(redshift)
    neg_redshift = - redshift[sort_indices]
    
    nper = n // nchunks
    indices = [0]
    for i in range(1, nchunks):
        isplit = np.searchsorted(neg_redshift, neg_redshift[indices[-1]+nper], side='right')
        indices.append(isplit)
    indices.append(n)
    print(f"Split indices: {indices}", flush=True)
    return indices
    

def init_bf(ofile, nchunks):
    bfw_list = []
    for i in range(nchunks):
        bf_w = BigFile(ofile + "-%d"%i, create=True)
        bfw_list.append(bf_w)
    return bfw_list


    
def io_col_data(bf_r, bfw_list, col, split_indices, sort_indices, nfiles):
    data_all = bf_r[col][:][sort_indices]
    for ichunk in range(nchunks):
        bf_w = bfw_list[ichunk]
        istart, iend = split_indices[ichunk], split_indices[ichunk+1]
        data = data_all[istart : iend]
        block = bf_w.create(col, name2type[col], len(data), nfiles)
        block.write(0, data)
        print(f"Finsihsed writing {col} to chunk {ichunk}", flush=True)


if __name__ == "__main__":

    # -------------- Cmd line Args ------------------------------
    parser = argparse.ArgumentParser(description="raw bhdetail to bigfile")
    parser.add_argument(
        "--ifile",
        required=True,
        type=str,
        help="path of the input file (end without /)",
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
    parser.add_argument(
        "--nchunks", required=True, type=int, help="chunks to split into"
    )

    args = parser.parse_args()
    nfiles = int(args.nfiles)
    nchunks = int(args.nchunks)
    bf_r = BigFile(args.ifile)

    redshift = bf_r['z'][:]
    sort_indices = get_sort_indices(redshift)
    print("got sorting indices", flush=True)
    split_indices = get_split_indices(redshift, sort_indices, nchunks)
    bfw_list = init_bf(args.ofile, nchunks)

    for col in SEL_COL:
        io_col_data(bf_r, bfw_list, col, split_indices, sort_indices, nfiles)
    print("Finished writing all columns", flush=True)
    
