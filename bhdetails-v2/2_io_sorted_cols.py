import numpy as np
from bigfile import BigFile
import struct
import glob
import os, sys
import argparse
import subprocess
import time


sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from sort_and_count import BlackholeDetails_BF
root = "/hildafs/datasets/Asterix/BHdetails-sorted"
sortidx_dir = root + "/sortidx"


def du(path):
    """disk usage in human readable format (e.g. '2,1GB')"""
    return subprocess.check_output(['du','-sh', path]).split()[0].decode('utf-8')


def io_col(file, snap):
    bf = BlackholeDetails_BF(file, snap)
    all_fields = bf.list_available_fields()
    print("Total number of cols: ", len(all_fields), flush=True)
    print("Total Size of this BF: ", bf.total_size_byte/1e9, "GBs", flush=True)

    sort_idx = np.load(sortidx_dir + f"/sort_idx_bhid_z_{snap:03d}.npy")

    # initialize new BF
    ofile = root + f"/BH-Details-R{snap:03d}"
    new_bf = BigFile(ofile, create=True)


    for field in all_fields:
        size = bf.get_field_size_byte(field)
        print(f"Size of {field}: {size/1e9} GBs", flush=True)

        time0 = time.time()
        data = bf.bf[field][:]
        time1 = time.time()
        print(f"Time taken to read: {time1-time0}", flush=True)


        time0 = time.time()
        data = data[sort_idx]
        time1 = time.time()
        print(f"Time taken to sort: {time1-time0}", flush=True)

        time0 = time.time()
        new_bf.create_from_array(field, data, Nfile = 16)
        time1 = time.time()
        print(f"Time taken to write: {time1-time0}", flush=True)

    # verify size of old and new BF
    print("Size of old BF: ", du(file), flush=True)
    print("Size of new BF: ", du(ofile), flush=True)



if __name__ == "__main__":

    # -------------- Cmd line Args ------------------------------
    parser = argparse.ArgumentParser(description="sort all cols by ID in a bhdetails file")
    parser.add_argument(
        "--snap",
        required=True,
        type=int,
        help="snap number",
    )

    parser.add_argument(
        "--path",
        required=True,
        type=str,
        help="path to the directory containing the old bhdetails files",
    )

    args = parser.parse_args()
    snap = args.snap
    file = args.path + f"/BH-Details-R{snap:03d}"

    io_col(file, snap)
    print("Done", flush=True)



    




    

