import numpy as np
from bigfile import BigFile,FileMPI
import h5py
import sys,os
import glob
import argparse
from bf_util import *
from mpi4py import MPI


def get_subfind_chunk(subroot):
    subdir = subroot + '/chunk*'
    
    chunk_list    = []
    maxgroup_list = []
    for ff in sorted(glob.glob(subdir)):
        cname    = ff.split('/')[-1]
        chunk    = int(cname.split('.')[0][5:])
        maxgroup = int(cname.split('.')[1])
        
        chunk_list.append(chunk)
        maxgroup_list.append(maxgroup)
        
    sort  = np.argsort(chunk_list)
    chunk_list    = np.array(chunk_list)[sort]
    maxgroup_list = np.array(maxgroup_list)[sort]
    return chunk_list, maxgroup_list 

def group_chunk_dir(groupidx):
    chunk = np.nonzero(maxgroup_list-1>=groupidx)[0][0]
    return chunk

    
def write_subgroup_sfr(block, gstart, gend):
    if gstart < 10000:
        for gidx in range(gstart,gend):
            start = dest_r['FOFGroups/OffsetByType'][gidx][0]
            end = start + dest_r['FOFGroups/LengthByType'][gidx][0]
            # data
            GroupSFR = dest_r['0/StarFormationRate'][start : end]
            
            FirstSubs = dest_r['FOFGroups/GroupFirstSub'][gidx]
            Nsubs     =  dest_r['FOFGroups/GroupNsubs'][gidx]
            SubLen    = dest_r['SubGroups/SubhaloLenType'][FirstSubs:FirstSubs+Nsubs][:,0].astype(np.int64)
            SubOff    = dest_r['SubGroups/SubhaloOffsetType'][FirstSubs:FirstSubs+Nsubs][:,0].astype(np.int64)
            
            SubOff = SubOff - SubOff[0:1]
            
            starts, ends = SubOff, SubOff + SubLen
            nsubs = len(SubOff)

            SubSFR = np.array([np.sum(GroupSFR[starts[i] : ends[i]]) for i in range(nsubs)])
            
            block.write(FirstSubs, SubSFR)
    else:
        start = dest_r['FOFGroups/OffsetByType'][gstart][0]
        end = dest_r['FOFGroups/OffsetByType'][gend][0]
        # data
        GroupSFR = dest_r['0/StarFormationRate'][start : end]
        
        FirstSubs = dest_r['FOFGroups/GroupFirstSub'][gstart : gend]
        Nsubs = dest_r['FOFGroups/GroupNsubs'][gstart : gend]

        hassub = (FirstSubs > 0).nonzero()[0]
        
        if len(hassub) == 0:
            return
        
        FirstSub = FirstSubs[hassub][0]
        LastSub  = FirstSubs[hassub][-1] + Nsubs[hassub][-1]

        SubLen    = dest_r['SubGroups/SubhaloLenType'][FirstSub : LastSub][:,0].astype(np.int64)
        SubOff    = dest_r['SubGroups/SubhaloOffsetType'][FirstSub : LastSub][:,0].astype(np.int64)
        
        
        SubOff = SubOff - SubOff[0]
        starts, ends = SubOff, SubOff + SubLen
        nsubs = len(SubOff)

        SubSFR = np.array([np.sum(GroupSFR[starts[i] : ends[i]]) for i in range(nsubs)])
    
        block.write(FirstSub, SubSFR)
            

def init_sfrblock():
    dtype  = '<f4'
    dsize  = dest_r['SubGroups/SubhaloSFR'].size
    nfile = dest_r['SubGroups/SubhaloSFR'].Nfile

    blockname = 'SubGroups/SubhaloSFR'
    block = dest_w.create(blockname, dtype, dsize, nfile)

    if rank == 0:
        print('Initizlized block:', blockname,flush=True)
        print('size:',dsize, 'type:',dtype,flush=True)
    
    return block

#--------------------------------------------------------------------------
 
if __name__ == "__main__":
    #----------- MPI Init ----------------
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    #-------------- Cmd line Args ------------------------------
    parser = argparse.ArgumentParser(description='subgrpID-subfind')
    parser.add_argument('--dest',required=True,type=str,help='path of the output file directory')
    parser.add_argument('--subroot',required=True,type=str,help='path of the subfind directory')
    parser.add_argument('--minpart',required=True,type=int,help='min particle in halo to start rewriting')
    # parser.add_argument('--gend',required=True,type=int,help='min particle in halo to start rewriting')

    args = parser.parse_args()
    
    #--------------------------
    minpart = int(args.minpart)
    subroot = args.subroot
    dest_w  = FileMPI(comm, args.dest, create=True)
    dest_r  = BigFile(args.dest)
    gsize   = dest_r['FOFGroups/LengthByType'].size
    ssize   = dest_r['SubGroups/SubhaloLenType'].size
    FUZZ    = 100000000
    
    Length  = dest_r['FOFGroups/LengthByType']
    Offset  = dest_r['FOFGroups/OffsetByType']
    
    #-------- find gend ----------------------
    gend = 0
    if rank == 0:
        LengthByType = Length[:]
        N01 = (LengthByType[:,0] + LengthByType[:,1]).astype(int)
        print("Minimum 0 and 1 parts in all groups:", min(N01), flush=True)
        gend = np.searchsorted(-N01, -minpart)
        print('Total Number of Groups:', Length.size,flush=True)
        del LengthByType, N01
        print('Ending Gidx for reordering: %d'%gend,flush=True)
        
    comm.barrier()
        
    gend = comm.allreduce(gend)
    
    comm.barrier()
    
    # ----------- Split tasks --------------------
    chunk_list, maxgroup_list = get_subfind_chunk(subroot)
    if minpart <= 32:
        cend = len(chunk_list)
    else:
        cend = group_chunk_dir(gend) + 1
        gend = maxgroup_list[cend]
        if rank == 1:
            print('Updated ending Gidx for reordering: %d'%gend,flush=True)

    #---------- Initialize all blocks --------------
    block = init_sfrblock()
    comm.barrier()
    
    Ngroups = gend
    istart = Ngroups * rank // size
    iend = Ngroups * (rank+1) // size
    
    print('rank %d processing groups %d to %d'%(rank,istart,iend),flush=True)
    comm.barrier()
    


    write_subgroup_sfr(block, istart, iend)
    
    