import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from bigfile import BigFile
import glob,os,struct
import argparse
from mpi4py import MPI

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
    isnap = (allsnaps["snap_num"] == sthis).nonzero()[0][0]
    
    if isnap == 0: # this is the first snap with 4/5
        snext = allsnaps[1]["snap_num"]
        znext = allsnaps[1]["redshift"]
        sprev = 0
        zprev = 100.
    elif isnap == len(allsnaps) - 1: # this is the last snap with 4/5
        snext = 1024
        znext = 0.
        sprev = allsnaps[-2]["snap_num"]
        zprev = allsnaps[-2]["redshift"]
    else:
        snext = allsnaps[isnap + 1]["snap_num"]
        sprev = allsnaps[isnap - 1]["snap_num"]
        znext = allsnaps[isnap + 1]["redshift"]
        zprev = allsnaps[isnap - 1]["redshift"]
    return zprev, znext

        
def select_mergers(mergers, zmax, zmin):
    print('Selecting mergers to process...')
    mask = mergers['z'] >= zmin
    mask &= mergers['z'] < zmax
    merger = mergers[mask]
    return merger

def find_bidx(BHIDs, bhid):
    try:
        return np.where(BHIDs == bhid)[0][0]
    except:
        return -1

def get_bhidxlist(bf, bhidlist):
    BHIDs = bf['5/ID'][:]
    bhidxlist = [find_bidx(BHIDs, bhid) for bhid in bhidlist]
    return bhidxlist

def find_subidx(bf, bidxlist):
    soff5 = bf['Subhalo/SubhaloOffsetType'][:][:, 5]
    mask = (bidxlist != -1).nonzero()[0]
    sidxlist = np.zeros(len(bidxlist), dtype=int)
    sidxlist_valid = np.searchsorted(soff5,bidxlist[mask],side='right')-1
    sidxlist[mask] = sidxlist_valid
    sidxlist[~mask] = -1
    return sidxlist


def append_field(array, field, value):
    if field in array.dtype.names:
        array[field] = value
    else:
        array = np.lib.recfunctions.append_fields(array, field, value, usemask=False)
    return array



if __name__ == "__main__":

    #-------------- Cmd line Args ------------------------------
    parser = argparse.ArgumentParser(description='subgrpID-subfind')
    parser.add_argument('--subroot',required=True,type=str,help='path of the subfind file directory')
    parse.add_argument('--snap',required=True,type=int,help='snapshot number')
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


    ibegin = rank * len(all_snap) // size
    iend = (rank + 1) * len(all_snap) // size
    if rank == size - 1:
        iend = len(all_snap)
    print(f"Rank {rank}: Processing snapshots from {ibegin} to {iend}", flush=True)
    for i in range(ibegin, iend):
        subdir = all_subdir[i]
        bf = BigFile(subdir)
        snap = all_snap[i]["snap_num"]
        zthis = all_snap[i]["redshift"]
        zprev, znext = find_zrange(snap)
        print(f"Rank {rank}: Processing snap {snap} with redshift {zthis}", flush=True)
        print(f"Delta z: {zprev} -> {zthis} -> {znext}", flush=True)

        #---------------- Process merger remnant for this snap ----------------
        merger_before = select_mergers(mergers, zmax=zprev, zmin=zthis)
        if len(merger_before) == 0:
            continue
        print(f"Rank {rank}: Processing snap {snap} with {len(merger_before)} mergers before this snap", flush=True)

        idlist1 = merger_before['id1']
        idlist2 = merger_before['id2']
        idrem = np.maximum(merger_before['id1'], merger_before['id2'])
        bhidxlist_rem = get_bhidxlist(bf, idrem)
        sidxlist_rem = find_subidx(bf, bhidxlist_rem)

        merger_before = append_field(merger_before, 'bhidx_rem', bhidxlist_rem)
        merger_before = append_field(merger_before, 'subidx_rem', sidxlist_rem)
        merger_before = append_field(merger_before, 'snap_rem', np.ones(len(merger_before), dtype=int) * snap)
        merger_before = append_field(merger_before, 'z_rem', np.ones(len(merger_before), dtype=float) * zthis)

        np.save(f"{outdir}/merger_rem_{snap}.npy", merger_before)
        print(f"Rank {rank}: Saved merger_rem_{snap}.npy", flush=True)
    comm.Barrier()

    if rank == 0:
        print("Finished processing merger remnant")
        # gather all merger_before
        merger_rem_all = []
        for snap in all_snap["snap_num"]:
            merger_rem = np.load(f"{outdir}/merger_rem_{snap}.npy")
            merger_rem_all.append(merger_rem)
        merger_rem_all = np.concatenate(merger_rem_all)
        np.save(f"{outdir}/merger_rem_all.npy", merger_rem_all)
    comm.Barrier()
    #--------------------- Now load again ------------------------
    mergers = np.load(f"{outdir}/merger_rem_all.npy")

    for i in range(ibegin, iend):
        #---------------- Process merger progenitor for this snap ----------------
        merger_after = select_mergers(mergers, zmax=zthis, zmin=znext)
        if len(merger_after) == 0:
            continue
        print(f"Rank {rank}: Processing snap {snap} with {len(merger_after)} mergers after this snap", flush=True)
        idlist1 = merger_after['id1']
        idlist2 = merger_after['id2']
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

        np.save(f"{outdir}/merger_after_{snap}.npy", merger_after)
        print(f"Rank {rank}: Saved merger_after_{snap}.npy", flush=True)

    comm.Barrier()
    if rank == 0:
        print("Finished processing merger progenitor")
        # gather all merger prog
        merger_all = []
        for snap in all_snap["snap_num"]:
            merger = np.load(f"{outdir}/merger_after_{snap}.npy")
            merger_all.append(merger)
        merger_all = np.concatenate(merger_all)
        np.save(f"{outdir}/merger_subidx_appended.npy", merger_all)

    comm.Barrier()
    if rank == 0:
        print("Finished appending subidx to merger catalog")
        



