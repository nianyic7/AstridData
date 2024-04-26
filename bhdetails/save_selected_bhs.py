import numpy as np
from bigfile import BigFile
import struct
import glob
import os, sys
from mpi4py import MPI
from bigfile import FileMPI
import argparse

SEL_COL=["acBHMass",  "acMass",  "BHID",  "BHMass",  "BHpos",  "BHvel",  "Mdot",  "Swallowed",  "SwallowID",  "z"]


def fetch_bh_info(bflist, bhid):
    dnames = SEL_COL
    dtypes = ["d",  "d", "q", "d", "3d", "3d",  "d", "i", "q", "d"]
    dt = np.dtype({"names": dnames, "formats": dtypes})
    bhinfo = []
    for bf in bflist:
        bhids = bf['BHID'][:]
        idx = np.where(bhids == bhid)[0]
        if len(idx) == 0:
            continue
        bhthis = np.zeros(len(idx), dtype=dt)
        for col in SEL_COL:
            bhthis[col] = bf[col][:][idx]
        bhinfo.append(bhthis)
    if len(bhinfo) == 0:
        return None
    bhinfo = np.concatenate(bhinfo)
    return bhinfo


def get_bflist(idir, smin, smax):
    filelist = []
    for f in sorted(glob.glob(idir + "/BH*")):
        snap = int(f.split('R')[-1][:3])
        if snap < smin or snap > smax:
            continue
        filelist.append(BigFile(f))
    return filelist


if __name__ == "__main__":
    # -------------- Cmd line Args ------------------------------
    parser = argparse.ArgumentParser(description="raw bhdetail to bigfile")
    parser.add_argument(
        "--idir",
        required=True,
        type=str,
        help="path of the input file (end without /)",
    )
    parser.add_argument(
        "--odir",
        required=True,
        type=str,
        help="path of the output file directory (end without /)",
    )
    parser.add_argument("--bhids", nargs="+", type=int, help="BHIDs to save", default=[])
    parser.add_argument(
        "--smin", type=int, help="min snapshot number", default=0
    )
    parser.add_argument(
        "--smax", type=int, help="max snapshot number", default=999
    )
    args = parser.parse_args()
    # ----------- MPI Init ----------------
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    idir = args.idir
    odir = args.odir
    bhids = args.bhids
    nbhs = len(bhids)
    bflist = get_bflist(idir, args.smin, args.smax)

    if rank == 0:
        print(f"Total number of BHs: {nbhs}", flush=True)
    ibegin = nbhs * rank // size
    iend = nbhs * (rank + 1) // size

    for i in range(ibegin, iend):
        print("rank %d processing BH %d" % (rank, bhids[i]), flush=True)
        bhid = bhids[i]
        bhinfo = fetch_bh_info(bflist, bhid)
        if bhinfo is None:
            print(f"BH {bhid} not found", flush=True)
            continue
        ofile = f"{odir}/BH-{bhid}.npy"
        np.save(ofile, bhinfo)
        print(f"Saved {ofile}", flush=True)
