import numpy as np
from bigfile import BigFile,FileMPI
import sys,os
import glob
import argparse
import matplotlib
import matplotlib.pyplot as plt
from mpi4py import MPI

if __name__ == "__main__":
    #----------- MPI Init ----------------
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    #-------------- Cmd line Args ------------------------------
    parser = argparse.ArgumentParser(description='subgrpID-subfind')
    parser.add_argument('--dest',required=True,type=str,help='path of the output file directory')
    parser.add_argument('--gstart',default=0,type=int,help='where to begin the rewriting')

    args = parser.parse_args()
    
    #--------------------------
    dest_w  = FileMPI(comm, args.dest, create=True)
    dest_r  = BigFile(args.dest)
    gstart  = int(args.gstart)    
    
    comm.barrier()
    #---------- Initialize  blocks --------------
    blockname = '5/SubgroupIndex'
    dtype = 'i8'
    dsize = dest_r['5/ID'].size
    nfile = dest_r['5/ID'].Nfile
    if gstart == 0:
        block5 = dest_w.create(blockname, dtype, dsize, nfile)
    else:
        block5 = dest_w[blockname]
    comm.barrier()    
    
    
    sLen = dest_r["SubGroups/SubhaloLenType"][:]
    sOff = dest_r["SubGroups/SubhaloOffsetType"][:]
    sOffbh = sOff[:,5]
    bh_num = len(dest_r["5/ID"][:])
    bh_idx_list = np.arange(bh_num)
    subhalo_search_idx = np.searchsorted(sOffbh[1:], bh_idx_list, side='right')
    
    grp_nsubs =dest_r['FOFGroups/GroupNsubs'][:]
    grp_firstsubs =dest_r['FOFGroups/GroupFirstSub'][:]
    bh_hostgrp = dest_r["5/GroupID"][:]-1
    grp_num = len(grp_nsubs)
    
    istart = grp_num * rank // size
    iend   = grp_num * (rank + 1) // size    
    
    comm.barrier()
    print('Rank %d starts from Group %d to Group %d'%(rank,istart,iend),flush=True)
    comm.barrier()
    
    test = []
    for grp_idx in range(3):
        
        first_sub_idx = grp_firstsubs[grp_idx]
        last_sub_idx = grp_firstsubs[grp_idx] + grp_nsubs[grp_idx]
        #bh_mask = subhalo_search_idx == last_sub_idx
        bh_mask_ingroup = bh_hostgrp == grp_idx
        bh_num_ingroup = np.sum(bh_mask_ingroup)
        bh_num_insubgroup = np.sum(sLen[first_sub_idx : last_sub_idx, 5])
        bh_num_wandering = bh_num_ingroup - bh_num_insubgroup
        
        if bh_num_wandering > 0:
            #print(subhalo_search_idx_test[bh_mask_ingroup][-bh_num_wandering:])
            #subhalo_search_idx_test = subhalo_search_idx.astype(float).copy()
            # subhalo_search_idx[bh_mask_ingroup][-bh_num_wandering:] = list((np.ones(bh_num_wandering)[:] * -1))
            # print(bh_num_wandering)
            # print(subhalo_search_idx[bh_mask_ingroup][-bh_num_wandering:])
            
            test.append(np.concatenate((subhalo_search_idx[bh_mask_ingroup][:-bh_num_wandering], np.ones(bh_num_wandering)[:] * -1)).astype(np.int64))    
        