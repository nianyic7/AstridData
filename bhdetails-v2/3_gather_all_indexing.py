import numpy as np
from bigfile import BigFile
import struct
import glob
import os, sys
import argparse
import subprocess
import time
import functools
import pickle

class SearchTable:
    def __init__(self, root = "/hildafs/datasets/Asterix/BHdetails-sorted/searchidx", minsnap=None, maxsnap=None):
        self.root = root
        self.minsnap = minsnap
        self.maxsnap = maxsnap
        self.snaps = self.get_search_snaps()
        self.snap2idx = {snap: idx for idx, snap in enumerate(self.snaps)}



    def get_search_snaps(self):
        self.files = sorted(glob.glob(self.root + "/search_idx*.npy"))
        snaps = []
        for file in self.files:
            snap = int(file.split("_")[-1].split(".")[0])
            if self.minsnap is not None and snap < self.minsnap:
                continue
            if self.maxsnap is not None and snap > self.maxsnap:
                continue
            snaps.append(snap)
        print("Total number of snaps: ", len(snaps), flush=True)
        print("first snap: ", min(snaps), flush=True)
        print("last snap: ", max(snaps), flush=True)
        return sorted(snaps)

    @functools.cached_property
    def all_ids(self):
        ids = set()
        for s in self.snaps:
            file = self.root + f"/search_idx_{s:03d}.npy"
            arr = np.load(file)
            bhid = arr["bhid"]
            ids.update(bhid)
            print(f"Loaded ids from snap: {file.split('_')[-1].split('.')[0]}", flush=True)
            print("Total number of ids: ", len(ids), flush=True)
        ids = np.array(list(sorted(ids)))
        print("Size of all ids: ", ids.nbytes/1e9, "GBs")
        mask = ids < 0
        print("Number of negative IDs: ", sum(mask))
        return ids

    
    def initialize_search_tables(self):
        print("Initializing search tables", flush=True)

        dlen = len(self.all_ids)
        itemsize = len(self.snaps)
        dnames = [str(snap) for snap in self.snaps]
        dtypes = [np.int64]*itemsize
        dtype = list(zip(dnames, dtypes))
        self.Length = np.zeros(dlen, dtype=dtype)
        self.Offset = np.zeros(dlen, dtype=dtype)


    def update_search_tables(self, snap):
        file = self.root + f"/search_idx_{snap:03d}.npy"
        arr = np.load(file)
        bhid = arr["bhid"]
        # print("min ID from search table:", min(bhid))

        idx = np.searchsorted(self.all_ids, bhid)
        self.Length[str(snap)][idx] = arr["len"]
        self.Offset[str(snap)][idx] = arr["off"]
        
        print(f"Updated search table for snap: {snap}", flush=True)

        # test
        idtest = 204265951658
        idxtest = np.where(self.all_ids == idtest)[0][0]
        print(f"test ID: {idtest} -> idx: {idxtest}")
        print("Length:", self.Length[str(snap)][idxtest])
        print("Offset:", self.Offset[str(snap)][idxtest])

    
    def save_search_tables_as_binary(self):
        self.Length.tofile(self.root + "/Length")
        self.Offset.tofile(self.root + "/Offset")
        print("Saved search tables as binary files", flush=True)

    def save_id2idx(self):
        id2idx = {id: idx for idx, id in enumerate(self.all_ids)}
        with open(self.root + "/id2idx.pkl", "wb") as f:
            pickle.dump(id2idx, f)
        print("Saved id2idx as pickle file", flush=True)

    

def excute(minsnap=None, maxsnap=None):
    # excute
    st = SearchTable(minsnap=minsnap, maxsnap=maxsnap)
    st.initialize_search_tables()
    for snap in st.snaps:
        st.update_search_tables(snap)
    st.save_search_tables_as_binary()
    st.save_id2idx()
    return st


def test(st, testids = None):
    print("Testing search tables....", flush=True)
    Nsnaps = len(st.snaps)
    SearchTableColNames = [str(snap) for snap in st.snaps]
    SearchTableColTypes = [np.int64] * Nsnaps
    SearchTableDataType = np.dtype(list(zip(SearchTableColNames, SearchTableColTypes)))
    SearchTableElemSize = np.zeros(1, dtype=SearchTableDataType).nbytes
    print("Size of each element in bytes:", SearchTableElemSize)

    id2idx = pickle.load(open(st.root + "/id2idx.pkl", "rb"))

    if testids is None:
        testids = np.random.choice(list(id2idx.keys()), 100)

    
    for ii in testids:
        print("reading ID:", ii)
        time0 = time.time()
        idx = id2idx[ii]
        offset = np.fromfile(st.root + "/Offset", dtype=SearchTableDataType, count=1, offset = idx * SearchTableElemSize)[0]
        length = np.fromfile(st.root + "/Length", dtype=SearchTableDataType, count=1, offset = idx * SearchTableElemSize)[0]

        for snap in st.snaps:
            if length[str(snap)] == 0:
                continue
            try:
                bf = BigFile(f"/hildafs/datasets/Asterix/BHdetails-sorted/BH-Details-R{snap:03d}")
            except:
                continue
            ids = bf["BHID"][offset[str(snap)]:offset[str(snap)]+length[str(snap)]]
            
            ids_read = set(ids)
            assert set(ids) == set([ii]), "IDs do not match, expected %d, got %d"%(ii, list(ids_read)[0])

        time1 = time.time()
    print("Passed!", flush=True)
    print("Time taken to read %d IDs from %d snaps: %f"%(len(testids), Nsnaps, time1-time0), flush=True)


if __name__ == "__main__":
    st = excute(minsnap=None, maxsnap=None)
    st = SearchTable()
    test(st)





