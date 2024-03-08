import numpy as np
from bigfile import BigFile
import struct
import glob
import os, sys
from mpi4py import MPI
from bigfile import FileMPI
import argparse

SEL_COL=['acBHMass', 'BHMass', 'Density', 'KEflag',  'Mdyn', 'NumDM', 'V1sumDM',
'acMass', 'BHpos', 'Encounter', 'KineticFdbkEnergy','MgasEnc', 'Swallowed', 'V2sumDM',
'BHID', 'BHvel', 'Entropy', 'DFAccel','DragAccel','GravAccel', 'Mdot', 'Mtrack', 'SwallowID',  'z']


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




        
def io_binary_file(ifile, ofile):
    dt = BHTYPE.TypeAll
    dtype = BHTYPE.dtype_all
    #-----------------------------------------------------------
    data = np.zeros(1, dtype=BHTYPE.TypeAll)
    nsingle = data.nbytes
    results = np.fromfile(ifile, dtype=dt, count=-1)
    hastrouble = 0
    
    mask =  (results["size"] != nsingle - 8) | (results["size2"] != nsingle - 8)
    mask = mask.nonzero()[0]
    if len(mask > 0):
        hastrouble = 1
        istart = mask[0]
        ntrouble = sum(mask)
        print(f"starting of trouble {istart}, total troubles {ntrouble}")
        offset = istart * nsingle
        # search for restart point
        for newoff in range(offset + 4, offset + nsingle):
            size1 = np.fromfile(ifile, dtype='i', count=1, offset=newoff)
            size2 = np.fromfile(ifile, dtype='i', count=1, offset=newoff+nsingle-4)
            if (size1 == nsingle - 8) & (size2 == nsingle - 8):
                print(f"resetting offset from {offset} to {newoff}")
                break
        next_data = np.fromfile(ifile, dtype=dt, count=1, offset=newoff)
        assert next_data["size"] == nsingle - 8, f"mismatch of size1 in file {ifile}, {next_data}"
        assert next_data["size2"] == nsingle - 8, f"mismatch of size2 in file {ifile}, {next_data}"
        results = np.fromfile(ifile, dtype=dt, count=-1, offset=newoff)
    
    results.tofile(ofile)






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
    # parser.add_argument(
    #     "--blocknames", nargs="+", help="blocks to save", default=None
    # )

    args = parser.parse_args()
    # ----------- MPI Init ----------------
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # --------------------------
    BHTYPE = BHType(SEL_COL)
    filelist = get_binary_files(args.ifile)
    Nraw = len(filelist)
    if rank == 0:
        print(f"Total number of files: {Nraw}", flush=True)
    comm.barrier()

    # split task and io
    fstart = Nraw * rank // size
    fend = Nraw * (rank + 1) // size
    print("Rank %03d will process file %03d to file %03d" % (rank, fstart, fend), flush=True)
    for i in range(fstart, fend):
        ifile = filelist[i]
        fname = ifile.split("/")[-1]
        ofile = args.ofile + "/" + fname
        io_binary_file(ifile, ofile)
    comm.barrier()
    if rank == 0:
        print("Done!", flush=True)