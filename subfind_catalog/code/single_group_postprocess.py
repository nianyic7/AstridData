import numpy as np
from bigfile import BigFile,FileMPI
import h5py as h5
import sys,os
import glob
import shutil
import argparse
from collections.abc import Iterable
import matplotlib
import matplotlib.pyplot as plt


def plot_subhalo(ax, itar, tab, grp, sOff, sLen, box=250000, zerocenter=True):
    beg, end = sOff[itar], sOff[itar] + sLen[itar]
    
    minpos = tab['Subhalo']['SubhaloPos'][itar]
    if zerocenter:
        pos1 = grp['PartType1']['Coordinates'][beg[1] : end[1]] - minpos
        pos4 = grp['PartType4']['Coordinates'][beg[4] : end[4]] - minpos
        pos1[pos1 > box/2] -= box
        pos1[pos1 < -box/2] += box
        pos4[pos4 > box/2] -= box
        pos4[pos4 < -box/2] += box
    else:
        pos1 = grp['PartType1']['Coordinates'][beg[1] : end[1]]
        pos4 = grp['PartType4']['Coordinates'][beg[4] : end[4]]
        
    
    sMass = tab['Subhalo']['SubhaloMassType'][:][:,4]
    mstar = sMass[itar]*1e10/0.6774
    group = tab['Subhalo']['SubhaloGroupNr'][:][itar]
    
    ax.scatter(pos1[:, 0], pos1[:, 1], label='DM', s=1)
    ax.scatter(pos4[:, 0], pos4[:, 1], label='Star', s=1)
    ax.scatter(0, 0, label='MinPos', color='cyan', s=50)
    #ax.set_title('Chunk%d, Sub%d, Mstar %.1e'%(chunk_idx, itar, mstar))    
    ax.set_title('Group %d Sub%d, Mstar %.1e'%(group, itar, mstar))  


def plot_subhalo_test(ax, itar, tab, grp, beg, end, box=250000, zerocenter=True):
    minpos = tab['Subhalo']['SubhaloPos'][itar]
    if zerocenter:
        pos1 = grp['PartType1']['Coordinates'][beg[1] : end[1]] - minpos
        pos4 = grp['PartType4']['Coordinates'][beg[4] : end[4]] - minpos
        pos1[pos1 > box/2] -= box
        pos1[pos1 < -box/2] += box
        pos4[pos4 > box/2] -= box
        pos4[pos4 < -box/2] += box
    else:
        pos1 = grp['PartType1']['Coordinates'][beg[1] : end[1]]
        pos4 = grp['PartType4']['Coordinates'][beg[4] : end[4]]


    pos1 = grp['PartType1']['Coordinates'][beg[1] : end[1]]
    pos4 = grp['PartType4']['Coordinates'][beg[4] : end[4]]        
    
    sMass = tab['Subhalo']['SubhaloMassType'][:][:,4]
    mstar = sMass[itar]*1e10/0.6774
    group = tab['Subhalo']['SubhaloGroupNr'][:][itar]
    
    ax.scatter(pos1[:, 0], pos1[:, 1], label='DM', s=1)
    ax.scatter(pos4[:, 0], pos4[:, 1], label='Star', s=1)
    ax.scatter(0, 0, label='MinPos', color='cyan', s=50)
    #ax.set_title('Chunk%d, Sub%d, Mstar %.1e'%(chunk_idx, itar, mstar))    
    ax.set_title('Group %d Sub%d, Mstar %.1e'%(group, itar, mstar))  
    
    
def write_tab_header(combine_tab, total_tab):
    parameter = combine_tab.create_group("Parameters")
    parameter_attrs = parameter.attrs
    for key in total_tab['Parameters'].attrs.keys():
        parameter_attrs[key] = total_tab['Parameters'].attrs[key]

    config = combine_tab.create_group("Config")
    config_attrs = config.attrs
    for key in total_tab['Config'].attrs.keys():
        config_attrs[key] = total_tab['Config'].attrs[key]

    ID_group = combine_tab.create_group("IDs")
    
    
    header = combine_tab.create_group("Header")
    header_attrs = header.attrs
    for key in total_tab['Header'].attrs.keys():
        if "Nsubhalos" in key:  # add separately
            continue 
        header_attrs[key] = total_tab['Header'].attrs[key]    
        
        
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
        
        
        
def copy_group_info(combine_tab, total_tab):
    Group_key = list(total_tab['Group'].keys())
    combine_tab.create_group("Group")
    
    for key in Group_key: 
        if key == 'GroupFirstSub':  # Added in the add_sub_info
            continue
        if key == 'GroupNsubs':  # Added in the add_sub_info
            continue
        # if key == "GroupOffsetType":  # Added separately 
        #     if isinstance(group_off, Iterable):
        #         data = total_tab['Group'][key][:]
        #         data[group_idx] = group_off
        #         combine_tab['Group'][key] = data
        #     else:
        #         combine_tab['Group'][key] = total_tab['Group'][key][:]                
        #     continue
            
        combine_tab['Group'][key] = total_tab['Group'][key][:]  
        
        
        
def add_sub_info(combine_tab, total_tab, single_tab, group_idx, group_off):
    group_idx_list = total_tab['Subhalo']['SubhaloGroupNr'][:]
    old_subgroup_num = np.sum(group_idx_list == group_idx)
    sub_start = total_tab['Group']['GroupFirstSub'][group_idx]
    sub_end = total_tab['Group']['GroupFirstSub'][group_idx+1]
    sub_num = single_tab['Header'].attrs['Nsubhalos_Total']
    
    combine_tab.create_group("Subhalo")
    for key in list(total_tab['Subhalo'].keys()):
        before_subhalo_data = total_tab['Subhalo'][key][:sub_start]
        after_subhalo_data = total_tab['Subhalo'][key][sub_end:]


        if key == 'SubhaloGroupNr':
            single_subhalo_data = np.ones(sub_num).astype("int64") * group_idx
        elif key == 'SubhaloOffsetType':
            single_subhalo_data = group_off + single_tab['Subhalo'][key][:]
            
        else:
            single_subhalo_data = single_tab['Subhalo'][key][:]

        if single_tab['Subhalo'][key][:].ndim > 1:
            data = np.vstack((before_subhalo_data, single_subhalo_data, after_subhalo_data))
        else:
            data = np.hstack((before_subhalo_data, single_subhalo_data, after_subhalo_data))
        
        combine_tab['Subhalo'][key] = data 


    # add 'group/groupfirstsub' and 'group/groupnsubs'
    group_nsubs = total_tab['Group/GroupNsubs'][:]
    group_nsubs[group_idx] = sub_num
    combine_tab['Group/GroupNsubs'] = group_nsubs

    first_sub = np.concatenate(([0], group_nsubs.cumsum()), dtype='i8', axis=0)
    first_sub = first_sub[:-1]
    
    nosub_group_mask = group_nsubs == 0
    first_sub[nosub_group_mask] = -1
    combine_tab['Group/GroupFirstSub'] = first_sub              
    
    
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='subgrpID-subfind')
    parser.add_argument('--snap_idx',required=True, type=int,help='snapshot idx')
    parser.add_argument('--chunk_idx',required=True, type=int,help='chunk idx')
    parser.add_argument('--subgroup_idx',required=True, type=int,help='subhalo idx')
    parser.add_argument('--subfind_datadir',required=True,type=str,help='path of the subfind file directory')
    #parser.add_argument('--savedir',required=True,type=str,help='path of save the single-group snapshot')
    parser.add_argument('--partition',required=True,type=str,help='path of saved partition numpy array')
    args = parser.parse_args()  
    
    
    print(args, flush=True)
    
    snap_idx = int(args.snap_idx)
    chunk_idx = int(args.chunk_idx)
    subgroup_idx = int(args.subgroup_idx)
    subfind_datadir = args.subfind_datadir
    partition = args.partition
    #singlegroup_dir = args.savedir
    
    
    max_group_list = np.load(partition)
    max_group = max_group_list[chunk_idx+1]
    datadir = os.path.join(subfind_datadir, f"chunk{chunk_idx}.{max_group}")
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
    
    beg, end = gOff[group_idx], gOff[group_idx] + gLen[group_idx]
    old_smass = sMass[group_idx]


    fig, axes = plt.subplots(1, 2)
    ax = axes[0]
    beg, end = gOff[group_idx], gOff[group_idx] + gLen[group_idx]
    print(beg, end)
        
    minpos = gPos[group_idx]
    print(minpos)
    pos1 = grp['PartType1']['Coordinates'][beg[1] : end[1]] - minpos
    pos4 = grp['PartType4']['Coordinates'][beg[4] : end[4]] - minpos
    box = 250000        
    pos1[pos1 > box/2] -= box
    pos1[pos1 < -box/2] += box
    pos4[pos4 > box/2] -= box
    pos4[pos4 < -box/2] += box


    ax.scatter(pos1[:, 0], pos1[:, 1], label='DM', s=1)
    ax.scatter(pos4[:, 0], pos4[:, 1], label='Star', s=1)
    ax.scatter(0, 0, label='MinPos', color='cyan', s=50)    

        
    snap_dir = f"/scratch3/09475/yihaoz/Astrid-subfind/{snap_idx}/single_group_{snap_idx}"
    datadir = os.path.join(snap_dir, f"chunk{chunk_idx}_{group_idx}/output/")
    tabfile = os.path.join(datadir, f"fof_subhalo_tab_{snap_idx}.hdf5")
    grpfile = os.path.join(datadir, f"snap_{snap_idx}.hdf5")

    tab = h5.File(tabfile, 'r')
    grp = h5.File(grpfile, 'r')

    sLen = tab['Subhalo']['SubhaloLenType'][:]
    sOff = tab['Subhalo']['SubhaloOffsetType'][:]
    cenID = tab['Subhalo']['SubhaloIDMostbound'][:]
    sMass = tab['Subhalo']['SubhaloMassType'][:][:,4]
    submass = tab['Subhalo/SubhaloMassType'][:]
    print(f"star in the subgroup: {submass[:,4].sum()}, star in the old group {old_smass}", flush=True)

    gLen = tab['Group']['GroupLenType'][:]
    gOff = tab['Group']['GroupOffsetType'][:]
    gPos = tab['Group']['GroupPos'][:]    
    #if(len(gLen) > 1):
    #    print("Warning! There are more than one group")    
    assert len(gLen) == 1, f"Warning! There are more than one group: {len(gLen)}!"
    
        

    ax = axes[1]

    #beg, end = gOff[group_idx], gOff[group_idx] + gLen[group_idx]
        
    minpos = gPos[0]
    pos1 = grp['PartType1']['Coordinates'][:] - minpos
    pos4 = grp['PartType4']['Coordinates'][:] - minpos
        
    pos1[pos1 > box/2] -= box
    pos1[pos1 < -box/2] += box
    pos4[pos4 > box/2] -= box
    pos4[pos4 < -box/2] += box

    

    ax.scatter(pos1[:, 0], pos1[:, 1], label='DM', s=1)
    ax.scatter(pos4[:, 0], pos4[:, 1], label='Star', s=1)
    ax.scatter(0, 0, label='MinPos', color='cyan', s=50)    
    plt.savefig(os.path.join(datadir, "single_group.png"))            
    
    
    prob_sub = 12
    #prob_sub = 30

    row = 5
    col = 5
    fig, axes = plt.subplots(row, col, figsize = (12, 12))

    for i in range(row):
        for j in range(col):
            plot_idx = i * col + j
            ax = axes[i][j]
            itar = plot_idx - row*col//2 + prob_sub
            plot_subhalo(ax, itar, tab, grp, sOff, sLen, box=250000, zerocenter=True)

    plt.subplots_adjust(wspace=0.3, hspace=0.4, top=0.93, bottom=0.02, left=0.02, right=0.93)
    fig.suptitle(f'chunk {chunk_idx}', fontsize=20)
    plt.savefig(os.path.join(datadir, "subhalo_beforecombine.png"))
    
    
        
    single_datadir = f"/scratch3/09475/yihaoz/Astrid-subfind/{snap_idx}/single_group_{snap_idx}/chunk{chunk_idx}_{group_idx}/output"
    print("single_datadir: ", single_datadir)
    tab_single = h5.File(os.path.join(single_datadir, f"fof_subhalo_tab_{snap_idx}.hdf5"))


    ## Need to be careful when load the original tab_file, especially when there are more than one group in this chunk that needs to be fixed 
    total_datadir = f"/scratch3/09475/yihaoz/Astrid-subfind/{snap_idx}/subfind_{snap_idx}/chunk{chunk_idx}.{max_group}/"
    print("total_datadir: ", total_datadir)

    if os.path.exists(os.path.join(total_datadir, "output-bk")):
        bk_filelist = os.listdir(os.path.join(total_datadir, "output-bk"))
        max_group_flag = -1
        for filename in bk_filelist:
            if f"fof_subhalo_tab_{snap_idx}_" in filename:
                print("There are one combined tab in the bk directory: ", filename)
                assert filename != f"fof_subhalo_tab_{snap_idx}_{group_idx}.hdf5",\
                    "The file corresponds to this group, needs to figure out where the proper total-tab is"
                print("previous group_idx:", filename.split("_")[-1].split(".")[0])
                group_flag = int(filename.split("_")[-1].split(".")[0])
                if group_flag > max_group_flag:
                    max_group_flag = group_flag
    else:
        print("Creating the back-up directory: ", os.path.join(total_datadir, "output-bk"))
        os.mkdir(os.path.join(total_datadir, "output-bk"))
        max_group_flag = -1

    if max_group_flag >= 0:
        tab_total_filename = os.path.join(total_datadir, "output-bk", f"fof_subhalo_tab_{snap_idx}_{max_group_flag}.hdf5")
        print("Loading tab file that fixed for group ", max_group_flag, ": ", tab_total_filename)
        tab_total = h5.File(tab_total_filename, 'r') 
    else:
        #tab_total_filename = os.path.join(total_datadir, "output", f"fof_subhalo_tab_{snap_idx}_trouble_{group_idx}.hdf5")
        tab_total_filename = os.path.join(total_datadir, "output", f"fof_subhalo_tab_{snap_idx}.hdf5")
        print("Loading original tab file: ", tab_total_filename, flush=True)
        tab_total = h5.File(tab_total_filename, 'r')     
        
        
    group_idx_list = tab_total['Subhalo']['SubhaloGroupNr'][:]
    old_subgroup_num = np.sum(group_idx_list == group_idx)

    newsub_num = tab_single['Header'].attrs['Nsubhalos_Total']
    print("old/new sub_num", old_subgroup_num, newsub_num, flush=True)

    sub_start = tab_total['Group']['GroupFirstSub'][group_idx]
    sub_end = tab_total['Group']['GroupFirstSub'][group_idx+1]

    ori_beg = tab_total['Group/GroupOffsetType'][group_idx]
    h5_newfile =  f"/scratch3/09475/yihaoz/Astrid-subfind/{snap_idx}/single_group_{snap_idx}/chunk{chunk_idx}_{group_idx}/output/fof_subhalo_tab_{snap_idx}_combinev2.hdf5"
    print("output file: ", h5_newfile, flush =True)
    
    
    with h5.File(h5_newfile,'w') as combine_h5: 
        write_tab_header(combine_h5, tab_total)
        combine_h5['Header'].attrs['Nsubhalos_Total'] = round(tab_total['Header'].attrs['Nsubhalos_Total'] + tab_single['Header'].attrs['Nsubhalos_Total'] - old_subgroup_num)
        combine_h5['Header'].attrs['Nsubhalos_ThisFile'] = round(tab_total['Header'].attrs['Nsubhalos_ThisFile'] + tab_single['Header'].attrs['Nsubhalos_ThisFile'] - old_subgroup_num)
        copy_group_info(combine_h5, tab_total)
        add_sub_info(combine_h5, tab_total, tab_single, group_idx, group_off=ori_beg)
        
        
    # do some check 
    with h5.File(h5_newfile,'r') as combine_h5: 
        for key in tab_total.keys():
            print(combine_h5[key], tab_total[key])

            attrs_keys_list = list(tab_total[key].attrs.keys())
            assert len(attrs_keys_list) == len(list(combine_h5[key].attrs.keys())), f"total_tab/{key} has {len(attrs_keys_list)} attrs while single_tab has {len(list(combine_h5[key].attrs.keys()))}"
            for attr_key in attrs_keys_list:
                if not combine_h5[key].attrs[attr_key] == tab_total[key].attrs[attr_key]:
                    print(f"Different attr: {key}.attrs:{attr_key} total_tab:{tab_total[key].attrs[attr_key]} single_tab:{combine_h5[key].attrs[attr_key]}")


        for group_key in tab_total['Group'].keys():
            old_data = tab_total['Group'][group_key][:]
            new_data = combine_h5['Group'][group_key][:]
            assert old_data.shape == new_data.shape, f"Group {group_key} has different shape in new_tab ({new_data.shape}) with in old_tab ({old_data.shape})"
            if not np.all(old_data == new_data):
                print(f"Different Group info: {group_key}")

            old_type = old_data.dtype
            new_type = new_data.dtype
            assert old_type == new_type, f"Group {group_key} has different type in new_tab ({new_type}) with in old_tab ({old_type})"

        new_tot_sub = newsub_num - old_subgroup_num + tab_total['Header'].attrs['Nsubhalos_Total']
        for sub_key in tab_total['Subhalo'].keys():
            old_data = tab_total['Subhalo'][sub_key][:]
            new_data = combine_h5['Subhalo'][sub_key][:]
            assert len(new_data) == new_tot_sub, f"new subhalo {sub_key} has {len(new_data)} subhalo_info, but not {new_tot_sub}!"
            
            
            
    print("Start reordering....")            
    snap_dir = f"/scratch3/09475/yihaoz/Astrid-subfind/{snap_idx}/subfind_{snap_idx}"
    output_bk = os.path.join(snap_dir, f"chunk{chunk_idx}.{max_group}/output-bk/")

    if max_group_flag >= 0:
        ori_snap_filename = os.path.join(output_bk, f"snap_{snap_idx}_{max_group_flag}.hdf5")
        ori_tab_filename = os.path.join(output_bk, f"fof_subhalo_tab_{snap_idx}_{max_group_flag}.hdf5")    
    else:
        ori_snap_filename = os.path.join(output_bk, f"snap_{snap_idx}.hdf5")
        ori_tab_filename = os.path.join(output_bk, f"fof_subhalo_tab_{snap_idx}.hdf5")

    print("ori_snap_filename: ", ori_snap_filename)
    print("ori_tab_filename: ", ori_tab_filename)
    print()

    #assert os.path.exists(ori_tab_filename) # need to check carefully which tab is the original one 

    if not os.path.exists(ori_tab_filename):
        ori_tab_address = os.path.join(snap_dir, f"chunk{chunk_idx}.{max_group}/output", f"fof_subhalo_tab_{snap_idx}.hdf5")
        print(f"back up the ori tab {ori_tab_address} to {ori_tab_filename}")
        shutil.copyfile(ori_tab_address, ori_tab_filename)
        

    if not os.path.exists(ori_snap_filename):
        
        ori_address = os.path.join(snap_dir, f"chunk{chunk_idx}.{max_group}/output", f"snap_{snap_idx}.hdf5")
        print(f"back up the ori snapshot {ori_address} to {ori_snap_filename}", flush=True)
        shutil.copyfile(ori_address, ori_snap_filename)


    ori_snap = h5.File(ori_snap_filename, 'r')
    ori_tab = h5.File(ori_tab_filename, 'r')    
                    
    single_dir = f"/scratch3/09475/yihaoz/Astrid-subfind/{snap_idx}/single_group_{snap_idx}"
    single_datadir = os.path.join(single_dir, f"chunk{chunk_idx}_{group_idx}/output/")
    single_tabfile = os.path.join(single_datadir, f"fof_subhalo_tab_{snap_idx}.hdf5")
    single_grpfile = os.path.join(single_datadir, f"snap_{snap_idx}.hdf5")

    single_tab = h5.File(single_tabfile, 'r')
    single_grp = h5.File(single_grpfile, 'r')            
    
    gOff = ori_tab['Group']['GroupOffsetType'][:]
    gLen = ori_tab['Group']['GroupLenType'][:]
    beg, end = gOff[group_idx], gOff[group_idx] + gLen[group_idx]


    group_partnum = single_grp['Header'].attrs["NumPart_Total"]

    assert np.all(single_grp['Header'].attrs["NumPart_Total"] == (end - beg))
    reorder_grp = os.path.join(single_datadir, f"snap_{snap_idx}_reorder.hdf5")
    print("reordered snapshot:", reorder_grp)

    astrid_filename = f"/scratch3/06431/yueyingn/ASTRID-II/output/PIG_{snap_idx}"
    bf_headfile = BigFile(astrid_filename)
    print('Astrid snapshot loaded', astrid_filename)        

    npart_tot = ori_snap['Header'].attrs["NumPart_Total"]

    with h5.File(reorder_grp,'w') as reorder_h5: 
        write_hdf_header(bf_headfile, reorder_h5, npart_tot, nfiles=0, npart_file=npart_tot)
        for ptype in [0,1,4,5]:
            reorder_h5.create_group("PartType" + str(ptype))

        for ptype in [0,1,4,5]:
            for part_key in ori_snap[f'PartType{ptype}'].keys():
                group_data = single_grp[f'PartType{ptype}/{part_key}'][:]
                ori_data = (ori_snap[f'PartType{ptype}/{part_key}'][:]).copy()
                ori_data[beg[ptype]:end[ptype]] = group_data

                reorder_h5[f'PartType{ptype}/{part_key}'] = ori_data
                print("Writing "+ "PartType"+str(ptype)+"/"+part_key, ori_data.shape)
            
            
    print(f"{snap_idx}, {chunk_idx}, {group_idx} Done")
    
    
        
    # check the combined tab
    snap_dir = f"/scratch3/09475/yihaoz/Astrid-subfind/{snap_idx}/single_group_{snap_idx}"
    datadir = os.path.join(snap_dir, f"chunk{chunk_idx}_{group_idx}/output/")

    tab_combine_file = h5_newfile
    grp_combine_file = os.path.join(datadir, f"snap_{snap_idx}_reorder.hdf5")
    print("tab combined file", tab_combine_file)
    print("group combined file", grp_combine_file, flush=True)


    tab_combine = h5.File(h5_newfile,'r')
    grp_combine = h5.File(grp_combine_file, 'r')

    sLen_combine = tab_combine['Subhalo']['SubhaloLenType'][:]
    sOff_combine = tab_combine['Subhalo']['SubhaloOffsetType'][:]    
    
    
    prob_sub = 12 + sub_start

    row = 5
    col = 5
    fig, axes = plt.subplots(row, col, figsize = (12, 12))

    for i in range(row):
        for j in range(col):
            plot_idx = i * col + j
            ax = axes[i][j]
            itar = plot_idx - row*col//2 + prob_sub
            plot_subhalo(ax, itar, tab_combine, grp_combine, sOff_combine, sLen_combine, box=250000, zerocenter=True)

    plt.subplots_adjust(wspace=0.3, hspace=0.4, top=0.93, bottom=0.02, left=0.02, right=0.93)
    fig.suptitle(f'chunk {chunk_idx}', fontsize=20)
    plt.savefig(os.path.join(datadir, "subhalo_aftercombine.png"))    
    
    
    # copy file
    output_bk_dir = os.path.join(f"/scratch3/09475/yihaoz/Astrid-subfind/{snap_idx}/subfind_{snap_idx}/chunk{chunk_idx}.{max_group}/output-bk")
    # make sure the original snapshot and tab_file (that contains problematic subhalos) are backed up
    assert os.path.exists(os.path.join(output_bk_dir, f"fof_subhalo_tab_{snap_idx}.hdf5")) 
    assert os.path.exists(os.path.join(output_bk_dir, f"snap_{snap_idx}.hdf5"))    
            
                
    snap_bk_file = os.path.join(output_bk_dir, f"snap_{snap_idx}_{group_idx}.hdf5")
    tab_bk_file = os.path.join(output_bk_dir, f"fof_subhalo_tab_{snap_idx}_{group_idx}.hdf5")

    print(shutil.copyfile(grp_combine_file, snap_bk_file))
    print(shutil.copyfile(tab_combine_file, tab_bk_file), flush=True)
                
        
        
                