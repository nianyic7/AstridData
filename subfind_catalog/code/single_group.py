import numpy as np
from bigfile import BigFile,FileMPI
import h5py as h5
import sys,os
import glob
import argparse
import matplotlib
import matplotlib.pyplot as plt

def write_hdf_header(bf, hdf5, npart_tot, nfiles, npart_file):
    """Generate an HDF5 header from a bigfile PIG header."""
    head = hdf5.create_group("Header")
    hattr = head.attrs
    battr = bf["Header"].attrs
    #As a relic from Gadget-1, the total particle numbers
    #are written as two separate 32 bit integers.
    hattr["NumPart_Total"] = np.uint32(npart_tot % 2**32)
    hattr["NumPart_Total_HighWord"] = np.uint32(npart_tot // 2**32)
    hattr["NumPart_ThisFile"] = np.int32(npart_file)
    hattr["NumFilesPerSnapshot"] = np.int32(nfiles)
    hattr["Redshift"] = 1./battr["Time"] - 1
    hattr["MassTable"] = np.int32([0,0,0,0,0,0]) # we input mass instead
    hattr["Flag_Cooling"] = 1
    hattr["Flag_Sfr"] = 1
    hattr["Flag_Feedback"] = 1
    hattr["Flag_DoublePrecision"] = 0
    hattr["Flag_StellarAge"] = 1
    hattr["Flag_Metals"] = 9
    #Pass other keys through unchanged. We whitelist expected keys to avoid confusing Gadget.
    hdfats = ["Time", "BoxSize", "Omega0", "OmegaLambda", "HubbleParam", "OmegaBaryon"]
    for attr in hdfats:
        hattr[attr] = battr[attr]
        
        
def plot_subhalo(ax, itar, tab, grp, sOff, sLen, box=250000):
    beg, end = sOff[itar], sOff[itar] + sLen[itar]
    
    minpos = tab['Subhalo']['SubhaloPos'][itar]
    pos1 = grp['PartType1']['Coordinates'][beg[1] : end[1]] - minpos
    pos4 = grp['PartType4']['Coordinates'][beg[4] : end[4]] - minpos
    sMass = tab['Subhalo']['SubhaloMassType'][:][:,4]
    
    pos1[pos1 > box/2] -= box
    pos1[pos1 < -box/2] += box
    pos4[pos4 > box/2] -= box
    pos4[pos4 < -box/2] += box
    
    mstar = sMass[itar]*1e10/0.6774
    
    ax.scatter(pos1[:, 0], pos1[:, 1], label='DM', s=1)
    ax.scatter(pos4[:, 0], pos4[:, 1], label='Star', s=1)
    ax.scatter(0, 0, label='MinPos', color='cyan', s=50)
    #ax.set_title('Chunk%d, Sub%d, Mstar %.1e'%(chunk_idx, itar, mstar))    
    ax.set_title('Sub%d, Mstar %.1e'%(itar, mstar))          
    
    
    
def draw_single_group(chunk_h5, hdf5name, header_bf, part_start, part_end, fnum, nfiles = 1):
    npart = part_end - part_start 
    
    startpart = part_start + fnum * (npart//nfiles)
    endpart = startpart + (npart//nfiles)    
    
    if fnum == nfiles - 1:
        endpart = part_end    
    
    #Open the file
    if nfiles==1:
        hfile = hdf5name+".hdf5"
    else:
        hfile = hdf5name + "."+str(fnum)+".hdf5"
        
    if os.path.exists(hfile):
        raise IOError("Not over-writing existing file ",hfile)

    
    with h5.File(hfile,'w') as hdf5:  
        #single_group_header(chunk_h5, hdf5, npart, nfiles, endpart - startpart)
        write_hdf_header(header_bf, hdf5, npart, nfiles, endpart - startpart)
        for ptype in [0,1,4,5]:
            hdf5.create_group("PartType" + str(ptype))
            
            for block in chunk_h5['PartType'+str(ptype)].keys():
                bfdata = chunk_h5['PartType'+str(ptype)][block][startpart[ptype]:endpart[ptype]]
                hdf5["PartType"+str(ptype)][block] = bfdata
                print("Writing "+ "PartType"+str(ptype)+"/"+block, bfdata.shape)    
                
                
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='subgrpID-subfind')
    parser.add_argument('--snap_idx',required=True, type=int,help='snapshot idx')
    parser.add_argument('--chunk_idx',required=True, type=int,help='chunk idx')
    parser.add_argument('--subgroup_idx',required=True, type=int,help='subhalo idx')
    parser.add_argument('--subfind_datadir',required=True,type=str,help='path of the subfind file directory')
    parser.add_argument('--savedir',required=True,type=str,help='path of save the single-group snapshot')
    parser.add_argument('--partition',required=True,type=str,help='path of saved partition numpy array')
    args = parser.parse_args()                
    
    print(args, flush=True)
    
    snap_idx = int(args.snap_idx)
    chunk_idx = int(args.chunk_idx)
    subgroup_idx = int(args.subgroup_idx)
    subfind_datadir = args.subfind_datadir
    partition = args.partition
    singlegroup_dir = args.savedir
    
    
    #max_group_list = np.load(f"/scratch3/09475/yihaoz/Astrid-subfind/{snap_idx}/idx_{snap_idx}.npy")
    #datadir = f"/scratch3/09475/yihaoz/Astrid-subfind/{snap_idx}/subfind_{snap_idx}/chunk{chunk_idx}.{max_group_idx}/"
    #savedir = f"/scratch3/09475/yihaoz/Astrid-subfind/{snap_idx}/single_group_{snap_idx}/chunk{chunk_idx}_{group_idx}/output"
    
    max_group_list = np.load(partition)
    max_group_idx = max_group_list[chunk_idx+1]
    datadir = os.path.join(subfind_datadir, f"chunk{chunk_idx}.{max_group_idx}")
    tabfile = os.path.join(datadir, "output", f"fof_subhalo_tab_{snap_idx}.hdf5")
    grpfile = os.path.join(datadir, "output", f"snap_{snap_idx}.hdf5")

    tab = h5.File(tabfile, 'r')
    grp = h5.File(grpfile, 'r')    
    
    
        
    gLen = tab['Group']['GroupLenType'][:]
    gOff = tab['Group']['GroupOffsetType'][:]
    gPos = tab['Group']['GroupPos'][:]
    sMass = tab['Group']['GroupMassType'][:][:,4]
    group_idx_list = tab['Subhalo']['SubhaloGroupNr'][:]    
    group_idx = group_idx_list[subgroup_idx]
    print("Problematic subgroup, group: ", subgroup_idx, group_idx, flush=True)
    
    #savedir = f"/scratch3/09475/yihaoz/Astrid-subfind/{snap_idx}/single_group_{snap_idx}/chunk{chunk_idx}_{group_idx}/output"
    savedir = os.path.join(singlegroup_dir, f"chunk{chunk_idx}_{group_idx}/output")
    if not os.path.exists(savedir):
        os.makedirs(savedir,)
        print("Creating", savedir)
        
        
    fig, ax = plt.subplots(1, 1)
    box=250000
    beg, end = gOff[group_idx], gOff[group_idx] + gLen[group_idx]
        
    minpos = gPos[group_idx]
    print(minpos)
    pos1 = grp['PartType1']['Coordinates'][beg[1] : end[1]] - minpos
    pos4 = grp['PartType4']['Coordinates'][beg[4] : end[4]] - minpos
        
    pos1[pos1 > box/2] -= box
    pos1[pos1 < -box/2] += box
    pos4[pos4 > box/2] -= box
    pos4[pos4 < -box/2] += box

    mstar = sMass[group_idx]*1e10/0.6774

    ax.scatter(pos1[:, 0], pos1[:, 1], label='DM', s=1)
    ax.scatter(pos4[:, 0], pos4[:, 1], label='Star', s=1)
    ax.scatter(0, 0, label='MinPos', color='cyan', s=50)        
        
        
    plt.savefig(os.path.join(savedir, "total_group.png"))
    
    
    sLen = tab['Subhalo']['SubhaloLenType'][:]
    sOff = tab['Subhalo']['SubhaloOffsetType'][:]
    cenID = tab['Subhalo']['SubhaloIDMostbound'][:]
    sMass = tab['Subhalo']['SubhaloMassType'][:][:,4]    
    
    prob_sub = subgroup_idx

    fig, axes = plt.subplots(3, 3, figsize = (10, 10))

    for i in range(3):
        for j in range(3):
            plot_idx = i * 3 + j
            ax = axes[i][j]
            itar = plot_idx - 4 + prob_sub
            plot_subhalo(ax, itar, tab, grp, sOff, sLen, box=250000)

    plt.subplots_adjust(wspace=0.2, hspace=0.2, top=0.93, bottom=0.02, left=0.02, right=0.93)
    fig.suptitle(f'chunk {chunk_idx}', fontsize=20)    
    
    plt.savefig(os.path.join(savedir, "problematic_subhalo.png"))    
    
    
    suboff = tab['Group/GroupFirstSub'][:]
    suboff = suboff[suboff >= 0]
    trouble_idx_list = np.nonzero(np.diff(tab['Subhalo/SubhaloLen'][:]) > 0)
    if len(trouble_idx_list[0]) > 0:
        check_idx = []
        for tr_idx in trouble_idx_list[0]:
            if not tr_idx+1 in suboff:
                check_idx.append(tr_idx)
        if len(check_idx):
            #print(chunk_idx, trouble_idx_list, check_idx) 
            print("checking problematic subhalos in this chunk:", chunk_idx, check_idx)     
            assert subgroup_idx in check_idx, f"{subgroup_idx} is not a problematic subhalo ({check_idx})."
            

    fig, axes = plt.subplots(3, 3, figsize = (10, 10))

    for i in range(3):
        for j in range(3):
            plot_idx = i * 3 + j
            ax = axes[i][j]
            itar = tab["Group/GroupFirstSub"][group_idx+plot_idx]
            plot_subhalo(ax, itar, tab, grp, sOff, sLen, box=250000)

    plt.subplots_adjust(wspace=0.2, hspace=0.2, top=0.93, bottom=0.02, left=0.02, right=0.93)
    fig.suptitle(f'chunk {chunk_idx}', fontsize=20)     
    plt.savefig(os.path.join(savedir, "first_subhalo_check.png"))               
    
    
    beg, end = gOff[group_idx], gOff[group_idx] + gLen[group_idx]
    astrid_filename = f"/scratch3/06431/yueyingn/ASTRID-II/output/PIG_{snap_idx}"
    bf_headfile = BigFile(astrid_filename)
    print('Astrid snapshot loaded', astrid_filename, flush=True)
    
        
    draw_single_group(chunk_h5=grp, hdf5name=os.path.join(savedir,f"snap_{snap_idx}"), header_bf=bf_headfile, 
                    part_start=beg, part_end=end, fnum=0, nfiles = 1)    