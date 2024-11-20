import numpy as np
from bigfile import BigFile
import struct
import glob
import os, sys
from mpi4py import MPI
from bigfile import FileMPI
import argparse

#SEL_COL=['acBHMass', 'BHMass', 'Density', 'KEflag',  'Mdyn', 'NumDM', 'V1sumDM',
#'acMass', 'BHpos', 'Encounter', 'KineticFdbkEnergy','MgasEnc', 'Swallowed', 'V2sumDM',
#'BHID', 'BHvel', 'Entropy', 'DFAccel','DragAccel','GravAccel', 'Mdot', 'Mtrack', 'SwallowID',  'z']
SEL_COL=["acMom",  "CountProgs", "DFAccel", "DragAccel",  "CountProgs", "Fdbk", "GasVel", 
"GravAccel", "MinPos", "MinPot",  "srParticles",  "srVel", "srDensity", "srDisp", "timebin"]

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

        if name_sel is None:
            name_sel = self.name_all
        self.name_sel = name_sel
        self.name2type = {name: dtype for name, dtype in zip(self.name_all, self.dtype_all)}


    @property
    def TypeAll(self):
        np_type = np.dtype({'names':self.name_all, 'formats':self.dtype_all})
        return np_type
    @property
    def TypeSel(self):
        name_sel = [name for name in self.name_all if name in name_sel]
        dtype_sel = [dtype for name, dtype in zip(self.name_all, self.dtype_all) if name in name_sel]
        np_type = np.dtype({'names':name_sel, 'formats':dtype_sel})
        return np_type


def get_binary_files(ifile):
    filelist = []
    for f in glob.glob(ifile + "/*"):
        filelist.append(f)
    return filelist


def get_length_offset(flist):
    Length = []
    Offset = []
    offset = 0
    data = np.zeros(1, dtype=BHTYPE.TypeAll)
    nsingle = data.nbytes
    print(f"Size of single element: {nsingle} bytes")

    for ff in flist[:]:
        file_stats = os.stat(ff)
        Length.append(file_stats.st_size // nsingle)
        Offset.append(offset)
        offset += Length[-1]
    return Length, Offset


        
def read_binary_file(file):
    dt = BHTYPE.TypeAll
    dtype = BHTYPE.dtype_all
    #-----------------------------------------------------------
    struct_fmt = ''.join(list(dtype))
    struct_len = struct.calcsize(struct_fmt)
    results = np.fromfile(file, dtype=dt, count=-1)

    results['z'] = 1./results['z'] - 1
    print("rank %d: read %d BHs"%(rank, len(results)))
    assert len(results) == Length[i], "number of BHs read in is not consistent with pre-determined size!"
    print(f"Total memory occupied: {results.nbytes/1e9} GBs")
    return results


def init_block(bf_w, blocknames, nfile):
    block_dict = {}
    for blockname in blocknames:
        dtype = BHTYPE.name2type[blockname]
        dsize = NTOT
        block = bf_w.create(blockname, dtype, dsize, nfile)
        block_dict[blockname] = block

        if rank == 0:
            print("Initizlized block:", blockname, flush=True)
    return block_dict

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
    # parser.add_argument(
    #     "--blocknames", nargs="+", help="blocks to save", default=None
    # )

    args = parser.parse_args()
    # ----------- MPI Init ----------------
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    bf_w = FileMPI(comm, args.ofile, create=True)
    nfiles = int(args.nfiles)
    

    # --------------------------
    BHTYPE = BHType(SEL_COL)
    filelist = get_binary_files(args.ifile)
    Nraw = len(filelist)
    if rank == 0:
        print(f"Total number of files: {Nraw}", flush=True)

    # get Nperfile
    Length, Offset = get_length_offset(filelist)
    NTOT = sum(Length)
    comm.barrier()
    if rank == 0:
        print(f"Total number of data: {NTOT}", flush=True)
    # initialize bf
    block_dict = init_block(bf_w, blocknames=SEL_COL, nfile=nfiles)
    comm.barrier()

    # split task and io
    fstart = Nraw * rank // size
    fend = Nraw * (rank + 1) // size
    print("Rank %03d will process file %03d to file %03d" % (rank, fstart, fend), flush=True)
    for i in range(fstart, fend):
        file = filelist[i]
        data = read_binary_file(file)
        for blockname in SEL_COL:
            block = block_dict[blockname]
            block.write(Offset[i], data[blockname])
        print("Rank %03d finished processing file %03d" % (rank, i), flush=True)
    comm.barrier()
    if rank == 0:
        print("Done!", flush=True)
