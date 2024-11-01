import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from bigfile import BigFile
import glob,os,struct
import argparse
from mpi4py import MPI
import time
from tqdm import tqdm
from numpy.lib import recfunctions as rfn

"""
Nianyi Chen 10/19/2024
For a given raw merger catalog, this code append the nearest subfind snapshot before and after the merger.
It also append the subhalo idx of the two host galaxies before the merger, and the remnant host galaxy after the merger.
"""

def get_all_subfind(root):
    all_subdir = sorted(glob.glob(root + '/PIG_*_subfind'))

    snaplist = []
    zlist = []
    ikeep = []
    for i, ss in enumerate(all_subdir):
        snap = int(ss.split('_')[-2])
        if snap in [225, 252, 272]:
            continue
        bf = BigFile(ss)
        scalefac = bf['Header'].attrs['Time'][0]
        reshift = 1./scalefac - 1
        snaplist.append(snap)
        zlist.append(reshift)
        ikeep.append(i)
    all_snap = np.zeros(len(snaplist), dtype=[('snap_num', 'i4'), ('redshift', 'f4')])
    all_snap['snap_num'] = snaplist
    all_snap['redshift'] = zlist
    all_subdir = [all_subdir[i] for i in ikeep]

    return all_subdir, all_snap


def find_zrange(sthis):
    isnap = (all_snap["snap_num"] == sthis).nonzero()[0][0]
    
    if isnap == 0: # this is the first snap with 4/5
        snext = all_snap[1]["snap_num"]
        znext = all_snap[1]["redshift"]
        sprev = 0
        zprev = 100.
    elif isnap == len(all_snap) - 1: # this is the last snap with 4/5
        snext = 1024
        znext = 0.
        sprev = all_snap[-2]["snap_num"]
        zprev = all_snap[-2]["redshift"]
    else:
        snext = all_snap[isnap + 1]["snap_num"]
        sprev = all_snap[isnap - 1]["snap_num"]
        znext = all_snap[isnap + 1]["redshift"]
        zprev = all_snap[isnap - 1]["redshift"]
    return zprev, znext

        
def select_mergers(mergers, zmax, zmin):
    print('Selecting mergers to process...')
    mask = mergers['z'] >= zmin
    mask &= mergers['z'] < zmax
    merger = mergers[mask]
    return merger

def find_bidx(id2idx, bhid):
    try:
        return id2idx[bhid]
    except KeyError:
        return -1

def get_bhidxlist(bf, bhidlist):
    BHIDs = bf['5/ID'][:]
    id2idx = {BHIDs[i]: i for i in range(len(BHIDs))}
    bhidxlist = [find_bidx(id2idx, bhid) for bhid in tqdm(bhidlist, desc='get_bhidxlist')]
    print(f"Length of bhidxlist: {len(bhidxlist)}")
    return np.array(bhidxlist)

def find_subidx(bf, bidxlist):
    soff5 = bf['SubGroups/SubhaloOffsetType'][:][:, 5]
    mask = (bidxlist>=0).nonzero()[0]
    sidxlist = np.zeros(len(bidxlist), dtype=int)
    sidxlist_valid = np.searchsorted(soff5,bidxlist[mask],side='right')-1
    sidxlist[mask] = sidxlist_valid
    sidxlist[~mask] = -1
    return sidxlist


def append_field(array, field, value):
    if field in array.dtype.names:
        array[field] = value
    else:
        array = rfn.append_fields(array, field, value, usemask=False)
    return array



if __name__ == "__main__":

    #-------------- Cmd line Args ------------------------------
    parser = argparse.ArgumentParser(description='subgrpID-subfind')
    parser.add_argument('--subroot',required=True,type=str,help='path of the subfind file directory')
    parser.add_argument('--snap',required=True,type=int,help='snapshot number')
    parser.add_argument('--mfile',required=True,type=str,help='raw merger npy file')
    parser.add_argument('--outdir',required=True,type=str,help='path of the output file directory')
    args = parser.parse_args()
    #--------------------------
    subroot = args.subroot
    mfile = args.mfile
    outdir = args.outdir
    snap = args.snap

    all_subdir, all_snap = get_all_subfind(subroot)

    print(f"Total number of subfind snapshots: {len(all_snap)}", flush=True)
    mergers = np.load(mfile)
    print(f"Total number of raw mergers: {len(mergers)}", flush=True)

    i = all_snap["snap_num"].tolist().index(snap)
    subdir = all_subdir[i]
    bf = BigFile(subdir)
    snap = all_snap[i]["snap_num"]
    zthis = all_snap[i]["redshift"]
    zprev, znext = find_zrange(snap)
    print(f"Delta z: {zprev} -> {zthis} -> {znext}", flush=True)

    #---------------- Process merger remnant for this snap ----------------
    t0 = time.time()
    merger_before = select_mergers(mergers, zmax=zprev, zmin=zthis)
    if len(merger_before) == 0:
        print(f"No mergers before snap {snap}", flush=True)
        exit()
    print(f"Processing snap {snap} with {len(merger_before)} mergers before this snap", flush=True)

    idlist1 = merger_before['ID1']
    idlist2 = merger_before['ID2']
    idrem = np.maximum(merger_before['ID1'], merger_before['ID2'])
    bhidxlist_rem = get_bhidxlist(bf, idrem)
    print("Got bhidxlist_rem")
    sidxlist_rem = find_subidx(bf, bhidxlist_rem)

    merger_before = append_field(merger_before, 'bhidx_rem', bhidxlist_rem)
    merger_before = append_field(merger_before, 'subidx_rem', sidxlist_rem)
    merger_before = append_field(merger_before, 'snap_rem', np.ones(len(merger_before), dtype=int) * snap)
    merger_before = append_field(merger_before, 'z_rem', np.ones(len(merger_before), dtype=float) * zthis)

    np.save(f"{outdir}/merger_rem_{snap}.npy", merger_before)
    print(f"Saved merger_rem_{snap}.npy", flush=True)


    #---------------- Process merger progenitor for this snap ----------------
    merger_after = select_mergers(mergers, zmax=zthis, zmin=znext)
    if len(merger_after) == 0:
        print(f"No mergers after snap {snap}", flush=True)
        exit()
    print(f"Processing snap {snap} with {len(merger_after)} mergers after this snap", flush=True)
    idlist1 = merger_after['ID1']
    idlist2 = merger_after['ID2']
    bhidxlist1 = get_bhidxlist(bf, idlist1)
    bhidxlist2 = get_bhidxlist(bf, idlist2)
    sidxlist1 = find_subidx(bf, bhidxlist1)
    sidxlist2 = find_subidx(bf, bhidxlist2)

    merger_after = append_field(merger_after, 'bhidx1', bhidxlist1)
    merger_after = append_field(merger_after, 'bhidx2', bhidxlist2)
    merger_after = append_field(merger_after, 'subidx1', sidxlist1)
    merger_after = append_field(merger_after, 'subidx2', sidxlist2)
    merger_after = append_field(merger_after, 'snap_prog', np.ones(len(merger_after), dtype=int) * snap)
    merger_after = append_field(merger_after, 'z_prog', np.ones(len(merger_after), dtype=float) * zthis)

    np.save(f"{outdir}/merger_prog_{snap}.npy", merger_after)
    print(f"Saved merger_prog_{snap}.npy", flush=True)

        



