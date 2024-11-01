import numpy as np
from bigfile import BigFile
import struct
import glob
import os, sys
import argparse
import time
import pickle

class BHinfo:
    def __init__(self, root="/hildafs/datasets/Asterix/BHdetails-sorted"):
        self.root = root
        self.idx_dir = root + "/searchidx"
        
        self.lenfile = self.idx_dir + "/Length"
        self.offfile = self.idx_dir + "/Offset"
        self.load_id2idx()
        self.snaps = self.get_search_snaps()
        self.snap2idx = {snap: idx for idx, snap in enumerate(self.snaps)}
        self.Nsnaps = len(self.snaps)


    def get_search_snaps(self):
        self.files = sorted(glob.glob(self.idx_dir + "/search_idx*.npy"))
        snaps = []
        for file in self.files:
            ss = int(file.split("_")[-1].split(".")[0])
            snaps.append(ss)
        print("Total number of snaps: ", len(snaps), flush=True)
        print("first snap: ", min(snaps), flush=True)
        print("last snap: ", max(snaps), flush=True)
        return sorted(snaps)
        
    def load_id2idx(self):
        with open(self.idx_dir + "/id2idx.pkl", 'rb') as f:
            self.id2idx = pickle.load(f)

    def load_bh_cols(self, bhidlist, cols):
        Nsnaps = len(self.snaps)
        SearchTableColNames = [str(snap) for snap in self.snaps]
        SearchTableColTypes = [np.int64] * Nsnaps
        SearchTableDataType = np.dtype(list(zip(SearchTableColNames, SearchTableColTypes)))
        SearchTableElemSize = np.zeros(1, dtype=SearchTableDataType).nbytes
        print("Size of each element in bytes:", SearchTableElemSize)

        # initialize data as dict for now
        all_data = dict()
        time0 = time.time()
        for ii in bhidlist:
            data = {c:[] for c in cols}
            print("reading ID:", ii)
            idx = self.id2idx[ii]
            print("index:", idx)
            length = np.fromfile(self.lenfile, dtype=SearchTableDataType, count=1, offset = idx * SearchTableElemSize)[0]
            offset = np.fromfile(self.offfile, dtype=SearchTableDataType, count=1, offset = idx * SearchTableElemSize)[0]
            for snap in self.snaps:
                if length[str(snap)] == 0:
                    continue
                try:
                    bf = BigFile(f"/hildafs/datasets/Asterix/BHdetails-sorted/BH-Details-R{snap:03d}")
                except Exception as e:
                    print(f"Cannot load bhdetails file BH-Details-R{snap:03d}: {e}")
                    sys.exit(1)
                ids = bf["BHID"][offset[str(snap)]:offset[str(snap)]+length[str(snap)]]
                ids_read = set(ids)
                assert set(ids) == set([ii]), "IDs do not match, expected %d, got %d"%(ii, list(ids_read)[0])
                
                ibeg = offset[str(snap)]
                iend = offset[str(snap)] + length[str(snap)]
                for c in cols:
                    data[c].append(bf[c][ibeg:iend])
            # turn dict into structured array
            # '3d' if shape is (N, 3), 'd' if dtype is float
            def set_dtype(c, data):
                if data[c].ndim == 1:
                    return data[c].dtype
                else:
                    return f"{data[c].shape[1]}d"
            for c in cols:
                data[c] = np.concatenate(data[c], axis=0)
            dtype = [(c, set_dtype(c, data)) for c in cols]
            data2 = np.array(list(zip(*[data[c] for c in cols])), dtype=dtype)
            all_data[ii] = data
        time1 = time.time()
        print(f"Time taken to read {len(bhidlist)} IDs from {Nsnaps} snaps: {time1-time0}", flush=True)
        return all_data

    def test(self, testids = None, cols = ["BHpos", "z", "BHID"]):
        # TODO: refine this test
        print("Testing...", flush=True)
        time0 = time.time()
        if testids is None:
            testids = np.random.choice(list(self.id2idx.keys()), 10)
        data = self.load_bh_cols(testids, cols)
        time1 = time.time()
        print("Passed!", flush=True)
        print("Time taken to read %d IDs: %f"%(len(testids), time1-time0), flush=True)

bhinfo = BHinfo()
bhinfo.test()


        