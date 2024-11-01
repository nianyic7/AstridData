import numpy as np
from bigfile import BigFile
import struct
import glob
import os, sys
import argparse
import functools
import time


"""
this script does the following:
1. read in the IDs of a bhdetail file, save as a set
"""

class BHType2:
    def __init__(self, name_sel=None):
        self.name_all = ('BHID','BHMass','Mdot','Density','timebin','Encounter','MinPos',\
             'MinPot','Entropy','GasVel','acMom','acMass','acBHMass',\
             'Fdbk','SPHID','SwallowID','CountProgs','Swallowed',\
             'BHpos','srDensity','srParticles','srVel','srDisp',\
             'DFAccel','DragAccel','GravAccel','BHvel','Mtrack','Mdyn',\
             'KineticFdbkEnergy','NumDM','V1sumDM','V2sumDM','MgasEnc','KEflag','z')
        self.dtype_all = ('q','d','d','d','i','i','3d',\
            'd','d','3d','3d','d','d',\
            'd','q','q','i','i',\
            '3d','d','d','3d','d',\
            '3d','3d','3d','3d','d','d',\
            'd','d','3d','d','d','i','d')

        if name_sel is None:
            name_sel = self.name_all
        self.name_sel = name_sel
        self.name2type = {name: dtype for name, dtype in zip(self.name_all, self.dtype_all)}
        self.TypeAll = np.dtype({'names':self.name_all, 'formats':self.dtype_all})


class BHType:
    def __init__(self, name_sel=None):
        self.name_all = ('BHID','BHMass','Mdot','Density','timebin','Encounter','MinPos',\
             'MinPot','Entropy','acMass','acBHMass',\
             'Fdbk','SPHID','SwallowID','CountProgs','Swallowed',\
             'BHpos','srDensity','srParticles','srVel','srDisp',\
             'DFAccel','DragAccel','GravAccel','BHvel','Mtrack','Mdyn','z')
        self.dtype_all = ('q','d','d','d','i','i','3d',\
            'd','d','3d','3d','d','d',\
            'd','q','q','i','i',\
            '3d','d','d','3d','d',\
            '3d','3d','3d','3d','d','d',\
            'd','d','3d','d','d','i','d')

        if name_sel is None:
            name_sel = self.name_all
        self.name_sel = name_sel
        self.name2type = {name: dtype for name, dtype in zip(self.name_all, self.dtype_all)}
        self.TypeAll = np.dtype({'names':self.name_all, 'formats':self.dtype_all})


class BlackholeDetails_BF:
    def __init__(self, path, snap):
        if snap >= 282:
            self.bhtype = BHType2()
        else:
            self.bhtype = BHType()
        self.bf = BigFile(path)
        self.BHIDs = None
        
    @property
    def data_length(self):
        return self.bf['BHID'].size

    @property
    def total_size_byte(self):
        length = self.data_length
        data = np.zeros(1, dtype=self.bhtype.TypeAll)
        nsingle = data.nbytes
        return length * nsingle

    def get_field_size_byte(self, field):
        field = self.bf[field]
        dsize = field.dtype.itemsize
        size_byte = dsize * field.size
        return size_byte

    def init_id_array(self):
        self.BHIDs = self.bf['BHID'][:]
        return

    def list_available_fields(self):
        return self.bf.keys()

    def get_data(self, field):
        return self.bf[field][:]

    @functools.cached_property
    def sort_idx_bhid_z(self):
        bhids = self.bf['BHID'][:]
        zz = self.bf['z'][:]
        data = np.column_stack((bhids, zz))
        sort_idx = np.lexsort((data[:, 1],data[:, 0]))
        return sort_idx

    def get_search_idx(self):
        """
        Group the data by BHID
        return: an array with bhid, ibeg, iend
        """
        sort_idx = self.sort_idx_bhid_z
        bhids = self.bf['BHID'][:]
        bhids = bhids[sort_idx]
        bhids = np.split(bhids, np.where(np.diff(bhids))[0] + 1)
        length = np.array([len(bhid) for bhid in bhids])
        unique_ids = np.array([bhid[0] for bhid in bhids])
        dtype = np.dtype([('bhid', 'i8'), ('off', 'i8'), ('len', 'i8')])
        data = np.zeros(length.size, dtype=dtype)
        data['bhid'] = unique_ids
        data['off'][1:] = np.cumsum(length)[:-1]
        data['len'] = length

        return data




if __name__ == "__main__":

    # -------------- Cmd line Args ------------------------------
    parser = argparse.ArgumentParser(description="save bh counters from each bhdetails file")
    parser.add_argument(
        "--sfirst",
        required=True,
        type=int,
        help="first snap number",
    )
    parser.add_argument(
        "--slast",
        required=True,
        type=int,
        help="last snap number",
    )
    parser.add_argument(
        "--path",
        required=True,
        type=str,
        help="path to the directory containing the bhdetails files",
    )

    odir = "/hildafs/datasets/Asterix/BHdetails-sorted"
    
    args = parser.parse_args()
    sfirst = args.sfirst
    slast = args.slast
    root = args.path


    snaplist = np.loadtxt("/hildafs/home/nianyic/AstridData/bhdetails-v2/all_file_nums.txt", dtype=int, delimiter=",")
    snaplist = snaplist[(snaplist >= sfirst) & (snaplist <= slast)]
    Nsnaptot = len(snaplist)
    print(f"Number of snapshots: {Nsnaptot}")

    for snap in snaplist:
        print(f"Processing snap: {snap}", flush=True)
        path = root + "/BH-Details-R%03d" % snap
        bh = BlackholeDetails_BF(path, snap)
        size = bh.total_size_byte
        print(f"Number of data elements: {bh.data_length}")

        time0 = time.time()
        sort_idx = bh.sort_idx_bhid_z
        time1 = time.time()
        print(f"Time taken to sort: {time1-time0}")
        print(f"Saving sorted index array to: {odir + '/sortidx/sort_idx_bhid_z_%03d.npy'%snap}")
        np.save(odir + "/sortidx/sort_idx_bhid_z_%03d.npy"%snap, sort_idx)

        time0 = time.time()
        search_idx = bh.get_search_idx()
        np.save(odir + "/searchidx/search_idx_%03d.npy"%snap, search_idx)
        time1 = time.time()
        print(f"Time taken to get search index: {time1-time0}")
        print(f"Saved: {odir + '/searchidx/search_idx_%03d.npy'%snap}")




    # print(f"Total size of the bigfile: {size/1e9} GBs")
    # size_pos = bh.get_field_size_byte('BHpos')
    # print(f"Size of BHPos: {size_pos/1e9} GBs")
    # all_fields = bh.list_available_fields()
    # print(all_fields)
    # print(len(all_fields))