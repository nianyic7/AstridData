from bigfile import BigFile
import numpy as np
import sys
import h5py
import sys,os
import pickle
import glob
import argparse

def get_obt(pig):
    """
    Args:
        pig (BigFile obj):               PIG file as BigFile
    Returns:
        OffsetByType (ndarray):          shape=(ngroup,6): [i,j] the starting index of species j in group i.
    """
    
    lbt = pig.open('FOFGroups/LengthByType')[:]
    obt = np.cumsum(lbt,axis=0)
    a1 = np.array([[0,0,0,0,0,0]],dtype=np.uint64)
    obt = np.append(a1,obt,axis=0)
    return obt

def get_avail_snaps(snapdir="/hildafs/datasets/Asterix/PIG*/PIG_*/"):
    """
    Args:
        snapdir (str):            generalized PIG file path [default: "/hildafs/datasets/Asterix/PIG*/PIG_*/"]
    Returns:
        snapz (ndarray):           (nsnap,2) array; snapz[:,0]: snapshot number, snapz[:,1]: snapshot redshift
    """
    
    
    pigfiles = sorted(list(set(glob.glob(snapdir))))
    snapz = []
    for d in pigfiles:
        pig = BigFile(d)
        scalefac = pig.open('Header').attrs['Time'][0]
        ii = int(d[-4:-1])
        zz = (1./scalefac - 1)
        snapz.append([ii,zz])
    
    return np.array(snapz)


def load_mergers(path,zmin=0,zmax=99):
    """
    Args:
        path (str):                  File(s) that stores the raw simulation mergers
        zmin (float):                min redshift cut [default: 0]
        zmax (float):                max redshift cut [default: 99]
    Returns:
        mergers (ndarray):           Structured array
    """
    
    binaries = []
    for f in sorted(glob.glob(path)):
        binaries.append(np.load(f))
    binaries = np.concatenate(binaries)
    print('Total number of mergers:',len(binaries),flush=True)
    
    mask = binaries['z'] >= zmin
    mask &= binaries['z'] <= zmax
    binaries = binaries[mask]
    print('Mergers after redshift cut:', len(binaries), flush=True)
    
    return binaries


def merger_to_snap(snapz,mergers):
    """
    Args:
        path (str):                  File(s) that stores the raw simulation mergers
        zmin (float):                min redshift cut [default: 0]
        zmax (float):                max redshift cut [default: 99]
    Returns:
        mergers (ndarray):           Structured array
    """
    
    binaries.sort(order='z')
    bins = np.digitize(binaries['z'],bins=snapz[1])
    mbinned = [binaries[bins==i] for i in range(len(snapz[1]))]
    blen = np.array([len(m) for m in mbinned])
    print('bin counts:',blen)
    
    
    

def ids_merged(mergers,snapz,i):
    """
    Find BHIDs of binaries that have merged just before this snapshot 
    snapz[i-1] > zmerge >= snapz[i]
    Since the merger already happened, we only save the larger ID

    Args:
        mergers (ndarray):        File(s) that stores the raw simulation mergers
        snapz (ndarray):          (2,nsnaps) array of all available snaps and their redshifts
        i (int):                  index of the snapshot in the array of snaps
    Returns:
        ids (set):                set of bhids to that has just merged by this snapshot
    """
    zmin = snapz[i][1]
    if i==0:
        zmax = 99.
    else:
        zmax = snapz[i-1][1]
        
    mask = mergers['z'] < zmax
    mask &= mergers['z'] >= zmin
    
    m = mergers[mask]
    ids = set(np.maximum(m['ID1'],m['ID2']))
    return ids
    
    
    
def ids_tomerge(mergers,snapz,i):
    """
    Find BHIDs of binaries that will merge in the next snapshot
    i.e. this is the last snapshot before the merger
    snapz[i] > zmerge >= snapz[i+1]
    we save both IDs since the merger has not happened

    Args:
        mergers (ndarray):        File(s) that stores the raw simulation mergers
        snapz (ndarray):          (nsnaps,2) array of all available snaps and their redshifts
        i (int):                  index of the snapshot in the array of snaps
    Returns:
        ids (set):                set of bhids to be merged soon after this snapshot
    """
    zmax = snapz[i][1]
    if i==len(snapz - 1):
        # FIXME: need to be careful
        # if the last merger is too far
        # from the last snapshot
        zmin = 0.
    else:
        zmin = snapz[i+1][1]
        
    mask = mergers['z'] < zmax
    mask &= mergers['z'] >= zmin
    
    m = mergers[mask]
    ids = set(np.concatenate([m['ID1'],m['ID2']]))
    return ids
    


def get_bh_groups(idset,pig):
    """
    Find BHIDs of binaries that will merge in the next snapshot
    i.e. this is the last snapshot before the merger
    snapz[i] > zmerge >= snapz[i+1]
    we save both IDs since the merger has not happened

    Args:
        mergers (ndarray):        File(s) that stores the raw simulation mergers
        snaps (ndarray):          (nsnaps,2) array of all available snaps and their redshifts
        i (int):                  index of the snapshot in the array of snaps
    Returns:
        ids (set):                set of bhids to be merged soon after this snapshot
    """

    # get obt
    lbt = pig['FOFGroups/LengthByType'][:]
    obt = get_obt(pig)
    
    # Load BHs
    bhgroup = pig['5/GroupID'][:]
    bhid = pig['5/ID'][:]
    
    groups = np.unique(np.array([bhgroup[i] for i,b in enumerate(bhid) if b in idset]))
    return groups


def get_pig(snap):
    if snap > 294:
        pig = BigFile(astrid_path + 'PIG2/PIG_%03d'%snap)
    else:
        pig = BigFile(astrid_path + 'PIG_files/PIG_%03d'%snap)
    return pig




#----------------------------------------------------------------------------

if __name__ == "__main__":
    
    # GLOBAL INFO
    astrid_path = '/hildafs/datasets/Asterix/'
    
    pig_path    = astrid_path + 'PIG*/'
    snapz       = get_avail_snaps(snapdir=pig_path+'PIG*/')
    
    merger_path = astrid_path + 'BH_details_bigfile/ASTRID-merger-catalog-z2.npy'
    mergers    = load_mergers(path=merger_path,zmin=0,zmax=99)
    
    #------------------- Arguments -------------------
    parser = argparse.ArgumentParser(description='save pre/post merger groups')
    parser.add_argument('--snap',required=True,type=int,help='snapshot number')
    parser.add_argument('--savedir',required=True,type=str,help='directory to save groupIDs')
    args = parser.parse_args()
    
    # ----------------- Do the Work ------------------
    snap     = int(args.snap)
    isnap    = (snapz[:,0]==snap).nonzero()[0][0]
    redshift = snapz[isnap][1]
    
    print('Processing Snapshot:',snap, 'Redshift:',redshift)
    pig      = get_pig(snap)

    idset  = ids_tomerge(mergers,snapz,isnap)
    idset.update(ids_merged(mergers,snapz,isnap))
    print('Unique BHs involved in mergers:',len(idset))

    groupIDs = get_bh_groups(idset,pig)

    print('Merger-relevant groups in snapshot %03d:'%snap, len(groupIDs),flush=True)
    
    np.save(args.savedir+'merger-group-%03d.npy'%snap, groupIDs)
    print('Done!',flush=True)



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

