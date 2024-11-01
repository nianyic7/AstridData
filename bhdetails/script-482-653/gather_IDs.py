import numpy as np
from bigfile import BigFile
import glob
import pickle

ROOT = "/hildafs/datasets/Asterix/BH_details_484-653/"
SMIN = 484
SMAX = 653
OUTDIR = "/hildafs/datasets/Asterix/BH_details_484-653"

def get_file_list():
    all_files = sorted(glob.glob(ROOT + "/BH-Details-R*.pkl"))
    flist = []
    for f in all_files:
        snum = int(f.split(".")[0].split("-")[-1][1:])
        if snum >= SMIN and snum <= SMAX:
            flist.append(f)
    return flist


def gather_IDs(flist):
    all_IDs = []
    for f in flist:
        print("processing:", f, flush=True)
        with open(f, "rb") as f:
            data = pickle.load(f)
        all_IDs += list(data.keys())
    all_IDs = np.array(list(set(all_IDs)))
    all_IDs = np.sort(all_IDs)
    return all_IDs

if __name__ == "__main__":
    flist = get_file_list()
    all_IDs = gather_IDs(flist)
    np.save(OUTDIR + "/ordered_all_IDs.npy", np.array(sorted(list(all_IDs))))
    print("Done!")
