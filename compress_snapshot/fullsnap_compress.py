import numpy as np
from bigfile import BigFile
from bigfile import FileMPI
import sys, os
from mpi4py import MPI
import argparse
import time
import glob

# we do not touch the integer columns

TypeMap = {
    np.dtype("float32"): np.dtype("float16"),
    np.dtype("float64"): np.dtype("float32"),
    np.dtype("uint64"): np.dtype("uint64"),
    np.dtype("uint32"): np.dtype("uint32"),
    np.dtype("uint8"): np.dtype("uint8"),
    np.dtype("int64"): np.dtype("int64"),
    np.dtype("int32"): np.dtype("int32"),
    np.dtype("int8"): np.dtype("int8"),
}

IntTypes = set(
    [
        np.dtype("uint64"),
        np.dtype("uint32"),
        np.dtype("uint8"),
        np.dtype("int64"),
        np.dtype("int32"),
        np.dtype("int8"),
    ]
)

# we also do not touch these cols as they may overflow
KeepCols = set(["ElectronAbundance", "Metallicity", "Metals"])


def compress(col, idata):
    """
    Transform and then compress data according to column properties
    and datatypes.
    Current transformation is 10th root then compress to ensure no overflow

    Args:
        col (str): data column name, e.g. "Position"
        idata (np.ndarray): array to be compressed

    Returns:
        np.ndarray: transformed array, same shape as idata

    """
    if idata.dtype in IntTypes:
        # leave integers alone
        return idata
    elif col == "Position":
        # direct compression into float32
        odata = idata.astype(TypeMap[idata.dtype])
        return odata
    elif col in KeepCols:
        # these columns may have overflow issue during compression
        # so leave them alone
        return idata
    else:
        # note: np.power(0.2) cannot deal with negative values
        dt = TypeMap[idata.dtype]
        odata = np.zeros_like(idata, dtype=dt)
        mask = idata > 0
        odata[mask] = np.power(idata[mask], 0.1).astype(dt)
        mask = idata < 0
        odata[mask] = -np.power(-idata[mask], 0.1).astype(dt)
        # note isfinite includes nan and inf checks
        mask_err = ~np.isfinite(odata)
        if mask_err.any():
            print(
                "Invalid data after compression in block %s, original data:" % col,
                idata[mask_err],
                flush=True,
            )

    return odata


def recover(col, odata):
    """
    Recover compressed data into original datatypes and values

    Args:
        col (str): data column name, e.g. "Position"
        odata (np.ndarray): array to be recovered

    Returns:
        np.ndarray: array recovered from compressed data
    """
    if odata.dtype in IntTypes:
        return odata
    elif col == "Position":
        rdata = odata.astype(np.dtype("float64"))
        return rdata
    elif col in KeepCols:
        return odata
    else:
        rdata = odata.astype(np.dtype("float32"))
        mask = rdata > 0
        rdata[mask] = np.power(rdata[mask], 10)
        mask = rdata < 0
        rdata[mask] = -np.power(rdata[mask], 10)

    return rdata


def rewrite_chunk(bf_r, bf_w, rank, Ntot, size, blockname, oblock):
    """
    compress and write a chunk of particle data

    Args:
        bf_r (bigfile): input bigfile
        bf_w (bigfile): output bigfile
        rank (int): this MPI rank
        Ntot (int): total length of this data block
        size (int): total MPI ranks
        blockname (str): particle type and data column, e.g. "1/Position"
        oblock (bigfile.datablock): bf col to write into, i.e. bf_w[blockname]
    """
    beg = Ntot * rank // size
    end = Ntot * (rank + 1) // size

    inner_size = nloops
    Nthis = end - beg
    col = blockname.split("/")[-1]

    for i in range(inner_size):
        ibeg = beg + Nthis * i // inner_size
        iend = beg + Nthis * (i + 1) // inner_size

        idata = bf_r[blockname][ibeg:iend]

        odata = compress(col, idata)
        oblock.write(ibeg, odata)

        # recover data and check for NaN and relative errors
        rdata = recover(col, odata)

        mask = idata != 0  # assume zeros are not messed up
        if not np.any(mask):
            continue
        # check NaN
        mask_err = ~np.isfinite(rdata)
        if mask_err.any():
            inan = idata[mask_err]
            onan = odata[mask_err]
            print(
                "WARNING: invalid data after recovery in col %s" % (blockname),
                flush=True,
            )
            print("Original data: ", inan[mask_err], flush=True)
            print("Compressed data:", onan[mask_err], flush=True)

        # check relative err within tolorance
        rel_err = np.abs(idata[mask] - rdata[mask]) / np.abs(idata[mask])
        maske = rel_err >= TOL
        if np.max(rel_err) >= TOL:
            print(
                "WARNING: error too large for %s: rel_error %.1e, original: %.1e, reconstruction: %.1e"
                % (
                    blockname,
                    np.max(rel_err),
                    idata[mask][maske][0],
                    rdata[mask][maske][0],
                ),
                flush=True,
            )
        # finish-up
        if i == inner_size:
            assert iend == end


def get_dtype_dsize(bf_r, blockname):
    """
    Get current and target datatype and data length
    of a block

    Args:
        bf_r (bigfile): input bigfile
        blockname (str): particle type and column name, e.g. "1/Position"

    Returns:
        dtype: _description_
        dtype:
        int :
    """
    old_type = bf_r[blockname][:0].dtype
    dshape = bf_r[blockname][:0].shape
    col = blockname.split("/")[-1]

    if len(dshape) == 1:  # simple dtype
        if col in KeepCols:
            new_type = old_type
        else:
            new_type = TypeMap[old_type]

    else:  # compound dtype
        dim = dshape[-1]
        if col in KeepCols:
            new_type = np.dtype((old_type, (dim,)))
        else:
            new_type = np.dtype((TypeMap[old_type], (dim,)))

    dsize = bf_r[blockname].size

    return old_type, new_type, dsize


def init_block(bf_r, bf_w, blockname, nfile):
    p = int(blockname.split("/")[0])
    old_type, dtype, dsize = get_dtype_dsize(bf_r, blockname)

    block = bf_w.create(blockname, dtype, dsize, nfile)

    if rank == 0:
        print("Initizlized block:", blockname, flush=True)
        print(
            "Nfile = {: <8d}, Ntot = {: <15d}, Old dtype = {: <10}".format(
                nfile, dsize, old_type.name
            ),
            "New dtype = ",
            dtype,
        )

    return block, dsize


if __name__ == "__main__":
    # ----------- MPI Init ----------------
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # -------------- Cmd line Args ------------------------------
    parser = argparse.ArgumentParser(description="PART file compression")
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
        "--tol", required=True, type=float, help="maximum relative error"
    )
    parser.add_argument(
        "--blocknames", nargs="+", help="blocks to process", default=None
    )

    args = parser.parse_args()

    # --------------------------

    bf_r = BigFile(args.ifile)
    bf_w = FileMPI(comm, args.ofile, create=True)
    nfiles = int(args.nfiles)
    nloops = int(args.nloops)
    TOL = args.tol
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
