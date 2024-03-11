import numpy as np
from bigfile import BigFile
import glob, os, struct
import warnings
import argparse
from find_single_pairs_astrid import *
from save_catalog_hdf5 import write_hdf5
from utils import *


if __name__ == "__main__":
    # ------------------- Arguments ------------------------------------
    parser = argparse.ArgumentParser(description="save-astrid-pair-cat")

    parser.add_argument("--snap", required=True, type=int, help="index of the PIG file")
    parser.add_argument("--cluster", required=True, type=str, help="cluster name")
    parser.add_argument("--odir", required=True, type=str, help="output directory")
    parser.add_argument(
        "--rmax", required=True, type=float, help="maximum separation for pair in kpc"
    )
    parser.add_argument(
        "--Mmin",
        default=0,
        type=float,
        help="minimum mass in Msun, no cut if not specified",
    )
    parser.add_argument(
        "--Lmin",
        default=0,
        type=float,
        help="minimum Lbol in erg/s, no cut if not specified",
    )

    args = parser.parse_args()
    snap = int(args.snap)
    rmax = args.rmax
    Mmin = args.Mmin
    Lmin = args.Lmin
    odir = args.odir
    DEBUG = True
    # ------------------ Read in BHs -----------------------------------
    if args.cluster == "vera":
        if snap <= 294:
            datadir = "/hildafs/datasets/Asterix/PIG_files"
        else:
            datadir = "/hildafs/datasets/Asterix/PIG2"

    elif args.cluster == "frontera":
        datadir = "/scratch3/06431/yueyingn/ASTRID-II/output"


    metadata = set_metadata_astrid()
    singles = load_singles_astrid(datadir, snap, Mmin, Lmin)
    pairs = find_pairs_astrid(singles, rmax)
    write_hdf5(f"{odir}/MBHpairs_{snap}.hdf5", metadata, singles, pairs, DEBUG=DEBUG)
    print(f"Done with snap {snap}!", flush=True)

    print("Done!", flush=True)
