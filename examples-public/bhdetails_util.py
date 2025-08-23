import numpy as np
from bigfile import BigFile
import struct
import glob
import os, sys
import argparse
import time
import pickle

'''
Note: This script is used to read BH details from the Astrid simulation.
For each BHdetails file, data are sorted by BHID then z.
Auxiliary files:
-  id2idx.pkl: dictionary mapping BHID to BH index in Length and Offset arrays
-  Length: (Nbhs, ) array, each element has datatype ["snap1", "snap2", ..., "Nsnap"]
    Length[i]["snap_j"] = number of entries for BHidx i in snap_i
-  Offset: (Nbhs, ) array, each element has datatype ["snap1", "snap2", ..., "Nsnap"]
    Offset[i]["snap_j"] = offset of BHidx i in snap_j
'''
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
        self.dtype_lean = ('i','q','f','f','f','i','i','3f',\
            'f','f','3f','3f','f','f',\
            'f','q','q','i','i',\
            '3f','f','f','3f','f',\
            '3f','3f','3f','3f','f','f',\
            'f','f','3f','f','f','i','f','i')
 
        if name_sel is None:
            name_sel = self.name_all
        self.name_sel = name_sel
        self.name2type = {name: dtype for name, dtype in zip(self.name_all, self.dtype_all)}
        self.name2type_lean = {name: dtype for name, dtype in zip(self.name_all, self.dtype_lean)}
 
 
    @property
    def TypeAll(self):
        np_type = np.dtype({'names':self.name_all, 'formats':self.dtype_all})
        return np_type
 
    @property
    def TypeLean(self):
        np_type = np.dtype({'names':self.name_all, 'formats':self.dtype_lean})
        return np_type
 
    @property
    def TypeSel(self):
        name_sel = [name for name in self.name_all if name in name_sel]
        dtype_sel = [dtype for name, dtype in zip(self.name_all, self.dtype_all) if name in name_sel]
        np_type = np.dtype({'names':name_sel, 'formats':dtype_sel})
        return np_type


#---------------------------------------------------------------------------------
class BHinfo:
    def __init__(self, root="/hildafs/datasets/Asterix/BHdetails-sorted"):
        self.root = root
        self.idx_dir = root + "/searchidx"
        
        self.lenfile = self.idx_dir + "/Length"
        self.offfile = self.idx_dir + "/Offset"
        self.load_id2idx()
        self.snaps = self.load_search_snaps()
        self.snap2idx = {snap: idx for idx, snap in enumerate(self.snaps)}
        self.Nsnaps = len(self.snaps)
        self.mergers = self.load_all_mergers()
        
        self.load_snap2z()
        self.BHType = BHType()

    def load_all_mergers(self, mfile=None):
        if mfile is None:
            mfile = "/hildafs/home/nianyic/scratch1/Astrid_data/lisa_mbhcat/up_to_date_mergers_nogal.npy"
        mergers = np.load(mfile)
        # sort from high redshift to low redshift
        mergers = mergers[np.argsort(-mergers["z"])]

        print("Total number of mergers: ", len(mergers), flush=True)
        zmin, zmax = np.min(mergers["z"]), np.max(mergers["z"])
        print("Redshift range: ", zmin, zmax, flush=True)
        return mergers

    def load_search_snaps(self):
        self.files = sorted(glob.glob(self.idx_dir + "/search_idx*.npy"))
        snaps = []
        for file in self.files:
            ss = int(file.split("_")[-1].split(".")[0])
            snaps.append(ss)
        print("Total number of snaps: ", len(snaps), flush=True)
        print("first snap: ", min(snaps), flush=True)
        print("last snap: ", max(snaps), flush=True)
        return sorted(snaps)

    def get_snap_zmin_zmax(self, snap):
        """
        Get the redshift range for a given snapshot
        and save it in self.snap2z if not already saved
        """
        if snap in self.snap2z:
            return self.snap2z[snap]
        bf = BigFile(f"/hildafs/datasets/Asterix/BHdetails-sorted/BH-Details-R{snap:03d}")
        z = bf["z"][:]
        zmin, zmax = np.min(z), np.max(z)
        self.snap2z[snap] = (zmin, zmax)
        return self.snap2z[snap]
    
    def load_snap2z(self):
        ofile = "/hildafs/datasets/Asterix/BHdetails-sorted/zrange.pkl"
        with open(ofile, "rb") as f:
            self.snap2z = pickle.load(f)
        
    
    def find_snapshot(self, zvalue):
        """
        Finds the snapshot number corresponding to a given z-value.

        Parameters:
            zvalue (float): The z-value to locate.
            snapshot_dict (dict): A dictionary mapping snapshot numbers to (zmin, zmax) tuples.

        Returns:
            int: The snapshot number containing the z-value.
            None: If the z-value is not in any range.
        """
        for snapshot, (zmin, zmax) in reversed(list(self.snap2z.items())):
            if zmin <= zvalue <= zmax:  # Check if zvalue is within the range
                return snapshot
        return None  # Return None if no match is found


        
    def load_id2idx(self):
        with open(self.idx_dir + "/id2idx.pkl", 'rb') as f:
            self.id2idx = pickle.load(f)


    def load_len_off_for_bhlist(self, bhidlist):
        Nsnaps = len(self.snaps)
        SearchTableColNames = [str(snap) for snap in self.snaps]
        SearchTableColTypes = [np.int64] * Nsnaps
        SearchTableDataType = np.dtype(list(zip(SearchTableColNames, SearchTableColTypes)))
        SearchTableElemSize = np.zeros(1, dtype=SearchTableDataType).nbytes
        print("Size of each element in bytes:", SearchTableElemSize)
        print("Loading Length and Offset for %d BHs..."%len(bhidlist), flush=True)
        Length = np.zeros((len(bhidlist),), dtype=SearchTableDataType)
        Offset = np.zeros((len(bhidlist),), dtype=SearchTableDataType)
        for i, bhid in enumerate(bhidlist):
            idx = self.id2idx[bhid]
            length = np.fromfile(self.lenfile, dtype=SearchTableDataType, count=1, offset = idx * SearchTableElemSize)[0]
            offset = np.fromfile(self.offfile, dtype=SearchTableDataType, count=1, offset = idx * SearchTableElemSize)[0]
            Length[i] = length
            Offset[i] = offset
        return Length, Offset


    def find_massive_prog_list(self, bhid, currlist):
        # Note: merger file already make m1 and id1 the massive BH
        # We are looking for cases where bhid is id2, and id2 is bigger than id1
        mask = (self.mergers["ID1"] == bhid) | (self.mergers["ID2"] == bhid)
        mergers = self.mergers[mask]
        print("Total mergers for bhid %d: %d"%(bhid, len(mergers)), flush=True)
        if len(mergers) == 0:
            # no merger at all
            return currlist
        mask = mergers["ID2"] == bhid
        mergers = mergers[mask]
        if len(mergers) == 0:
            # this bh is always the massive one
            print("bh %d is always the massive one"%bhid, flush=True)
            return currlist
        mask = mergers["ID2"] > mergers["ID1"]
        mergers = mergers[mask]
        if len(mergers) == 0:
            # this bh is swallowed by a bigger one already
            # should not happen if we are tracing an existing bh
            print("bh %d is swallowed by a bigger one already"%bhid, flush=True)
            return currlist
        most_recent = mergers[-1]
        idmassive = most_recent["ID1"]
        zmerge = most_recent["z"]
        currlist.append((idmassive, zmerge))
        result = self.find_massive_prog_list(idmassive, currlist)
        return result

    
    def get_bhids_massive_progs(self, bhidlist):
        """
        Get the IDs of all massive progs and current bhids for loading
        """
        massive_progs_dict = dict()
        for bhid in bhidlist:
            massive_progs_dict[bhid] = self.find_massive_prog_list(bhid, [])
        all_progs =[]
        for k, v in massive_progs_dict.items():
            progs = [x[0] for x in v]
            all_progs = all_progs + progs
        all_progs = np.unique(all_progs)
        bhs_to_load = np.unique(np.concatenate([bhidlist, all_progs]))
        return massive_progs_dict, bhs_to_load

    def find_prog_in_snap(self, bhid, snap, massive_progs_dict):
        """
        Find the most progenitor of bhid in a given snapshot
        Note: currently not used
        """
        if len(massive_progs_dict[bhid]) == 0:
            return bhid
        progs = massive_progs_dict[bhid]
        zmin, zmax = self.get_snap_zmin_zmax(snap)
        zmerge = [progs[i][1] for i in range(len(progs))]
        sort_idx = np.argsort(np.array(zmerge)) # redshift low to high
        pids = np.array([progs[i][0] for i in sort_idx])
        zmerge = np.array(zmerge)[sort_idx]

        flag = False
        progid_list = []
        for i, z in enumerate(zmerge):
            if not flag and (z >= zmin):
                flag = True
                if i == 0:
                    progid_list.append(bhid)
                else:
                    progid_list.append(pids[i-1])
            if flag:
                if z <= zmax:
                    progid_list.append(pids[i])
                else:
                    break
        
        return progid_list


    def load_bh_cols_massive_progs(self, bhidlist, cols, smin=0, smax=1024):
        """
        Load bhdetails not only for the given bhidlist, but also for all massive progenitors of these bhs.
        steps:
        1. recursively trace all mergers of the most massive progenitors, save as a dict
        2. load bhdetails for all bhs using load_bh_cols_simple
        """
        massive_progs_dict, bhs_to_load = self.get_bhids_massive_progs(bhidlist)
        all_data = self.load_bh_cols_simple(bhs_to_load, cols, smin, smax)
        return massive_progs_dict, all_data


    def load_bh_cols_simple(self, bhidlist, cols, smin=0, smax=1024):
        if "Swallowed" not in cols:
            cols.append("Swallowed")
        print("Total BHs to load:", len(bhidlist), flush=True)
        Length, Offset = self.load_len_off_for_bhlist(bhidlist)
        NtotBH = [sum(Length[i]) for i, bhid in enumerate(bhidlist)] # total data length of each bh
        time0 = time.time()
        
        
        all_data = {c:[] for c in cols}
        for snap in self.snaps:
            if (snap > smax) | (snap < smin):
                continue
            print("processing snap:", snap, flush=True)
            # if there is no data for any bhs in this snap, skip
            if np.sum(Length[str(snap)]) == 0:
                print("Skiping snap %d, no data for any BHs"%snap, flush=True)
                continue
            bf = BigFile(f"/hildafs/datasets/Asterix/BHdetails-sorted/BH-Details-R{snap:03d}")
            ibeg, iend = Offset[str(snap)], Offset[str(snap)] + Length[str(snap)]
            for c in cols:
                COL = bf[c]
                # load data for all bhs in one go
                data = [COL[ibeg[i] : iend[i]] for i, bhid in enumerate(bhidlist)]
                all_data[c].append(data)
                
            print("time taken: ", time.time() - time0, flush=True)
            time0 = time.time()
        
        for c in cols:
            lists = all_data[c] # a list of Nsnap lists, where the isnapth list has Nbh arrays
            all_data[c] = [np.concatenate(elements, axis=0) for elements in zip(*lists)] # a list of Nbh arrays
        
        dtype = [(c, self.BHType.name2type_lean[c]) for c in cols]
        bhdata_all = {bhid: np.zeros(NtotBH[i], dtype=dtype) for i, bhid in enumerate(bhidlist)}
        for i, bhid in enumerate(bhidlist):
            for c in cols:
                bhdata_all[bhid][c] = all_data[c][i]
            bhdata = bhdata_all[bhid]
            mask = bhdata["Swallowed"] == 0
            bhdata = bhdata[mask]
            sortidx = np.argsort(-bhdata["z"])
            bhdata = bhdata[sortidx]
            bhdata_all[bhid] = bhdata
            
        print("Finished loading data for %d BHs"%len(bhidlist), flush=True)
        t1 = time.time()
        print("Time taken: ", t1-time0, flush=True)
        return bhdata_all

    def test(self, testids = None, cols = ["BHpos", "z", "BHID"]):
        # TODO: refine this test
        print("Testing...", flush=True)
        time0 = time.time()
        if testids is None:
            testids = np.random.choice(list(self.id2idx.keys()), 10)
        # print("test simple loading...", flush=True)
        # data = self.load_bh_cols_simple(testids, cols)
        time1 = time.time()
        # print("Passed!", flush=True)
        # print("Time taken to read %d IDs: %f"%(len(testids), time1-time0), flush=True)

        print("test massive progs loading...", flush=True)
        prog_dict, data = self.load_bh_cols_massive_progs(testids, cols)
        for k, v in prog_dict.items():
            print(k, v)
        time2 = time.time()
        print("Passed!", flush=True)
        print("Time taken to read %d IDs: %f"%(len(testids), time2-time1), flush=True)



if __name__ == "__main__":
    bhinfo = BHinfo()
    bhinfo.test()


        
