import numpy as np
from bigfile import BigFile
import struct
import glob
import os, sys
from mpi4py import MPI




def get_binary_files(ifile):
    filelist = []
    for f in glob.glob(ifile + "/*"):
        filelist.append(f)
    print("Found %d binary files" % len(filelist))
    return filelist


        
def read_binary_file(file):
    fields = ('size','BHID','BHMass','Mdot','Density','timebin','encounter','MinPos',\
             'MinPot','Entropy','GasVel','acMom','acMass','acBHMass',\
             'Fdbk','SPHID','SwallowID','CountProgs','Swallowed',\
             'BHpos','srDensity','srParticles','srVel','srDisp',\
             'DFAccel','DragAccel','GravAccel','BHvel','Mtrack','Mdyn',\
             'KinFdbk','NumDM','V1sumDM','V2SumDM','MgasEnc','KEflag','z','size2')
    dtype = ('i','q','d','d','d','i','i','3d',\
            'd','d','3d','3d','d','d',\
            'd','q','q','i','i',\
            '3d','d','d','3d','d',\
            '3d','3d','3d','3d','d','d',\
            'd','d','3d','d','d','i','d','i')
    dt = {'names':fields, 'formats':dtype}
    #-----------------------------------------------------------
    struct_fmt = ''.join(list(dtype))
    struct_len = struct.calcsize(struct_fmt)
    results = []
    for filename in bhfile_zoo[istart:iend]:
        data = np.fromfile(filename, dtype=dt, count=-1)
        results.append(data)
    results = np.concatenate(results)

    results['z'] = 1./results['z'] - 1
    print("rank %d: read %d BHs from %d files"%(rank, len(results), len(bhfile_zoo[istart:iend])))
    print(f"Total memory occupied: {results.nbytes/1e9} GBs")
    return results












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
        "--nloops", required=True, type=int, help="number of loops per task"
    )
    parser.add_argument(
        "--blocknames", nargs="+", help="blocks to save", default=None
    )

    args = parser.parse_args()
    # ----------- MPI Init ----------------
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # --------------------------
    filelist = get_binary_files(args.ifile)
    Nraw = len(filelist)
    istart = Nraw * rank // size
    iend = Nraw * (rank + 1) // size
    print("Rank %03d will process file %03d to file %03d" % (rank, istart, iend), flush=True)



   








    bf_w = FileMPI(comm, args.ofile, create=True)
    nfiles = int(args.nfiles)
    nloops = int(args.nloops)

    # block names to process
    blocknames = []
    if args.blocknames is None:
        all_blocks = sorted(glob.glob(args.ifile + "/*/*"))
        for ff in all_blocks:
            column = ff.split("/")[-1]
            ptype = ff.split("/")[-2]
            blockname = ptype + "/" + column

            if ptype in ["0", "1", "4", "5"]:
                blocknames.append(blockname)
    else:
        for bn in args.blocknames:
            if os.path.isdir(args.ifile + "/" + bn):
                blocknames.append(bn)
            elif rank == 0:
                print("Block %s does not exit!" % (bn), flush=True)

    comm.barrier()

    if rank == 0:
        print("All blocks to process:", flush=True)
        for bb in blocknames:
            print(bb, flush=True)

    comm.barrier()
    if rank == 0:
        t1 = time.time()
    for blockname in blocknames:
        oblock, Ntot = init_block(bf_r, bf_w, blockname, nfiles)
        comm.barrier()
        rewrite_chunk(bf_r, bf_w, rank, Ntot, size, blockname, oblock)
        comm.barrier()

        if rank == 0:
            t2 = time.time()
            print("Finished Processing {: <15}".format(blockname), flush=True)
            print("time used = %.1f min:" % ((t2 - t1) / 60.0), flush=True)
            t1 = t2
        comm.barrier()

    if rank == 0:
        print("Done!", flush=True)
