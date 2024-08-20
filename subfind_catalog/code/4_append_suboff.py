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


def get_SubgroupOff(GroupOff,SubLen):
    if len(SubLen) == 0:
        return np.zeros((0,6),dtype = np.int64)
    SubOff = np.zeros_like(SubLen,dtype=np.int64)
    SubOff[1:] = np.cumsum(SubLen[:-1],axis=0)
    return (SubOff + GroupOff).astype(np.int64)

def write_subgroup_offset(block,gstart,gend):
    if gstart < 200:
        for gidx in range(gstart,gend):
            FirstSubs = dest_r['FOFGroups/GroupFirstSub'][gidx]
            Nsubs     =  dest_r['FOFGroups/GroupNsubs'][gidx]
            GroupOff  = dest_r['FOFGroups/OffsetByType'][gidx]
            SubLen    = dest_r['SubGroups/SubhaloLenType'][FirstSubs:FirstSubs+Nsubs]
            SubOff    = get_SubgroupOff(GroupOff,SubLen)
            block.write(FirstSubs,SubOff)
            # print('group %d done!'%gidx)
            
    else:
        FirstSubs = dest_r['FOFGroups/GroupFirstSub'][gstart:gend]
        mask = FirstSubs > 0
        FirstSubs = FirstSubs[mask]
        if len(FirstSubs) == 0:
            print(f"no subhalo in the group ({gstart} to {gend})!")
            return
        Nsubs     =  dest_r['FOFGroups/GroupNsubs'][gstart:gend][mask]
        GroupOff  = dest_r['FOFGroups/OffsetByType'][gstart:gend][mask]
    
    
        sstart,send = FirstSubs[0], FirstSubs[-1] + Nsubs[-1]
        try:
            SubLen    = dest_r['SubGroups/SubhaloLenType'][sstart:send]
        except ValueError:
            print(gstart,gend,sstart,send,flush=True)
            assert 0==1
        FirstSubs = FirstSubs - FirstSubs[0]

        starts, ends = FirstSubs, FirstSubs + Nsubs
        SubOff = [get_SubgroupOff(GroupOff[i],SubLen[starts[i]: ends[i]]) for i in range(len(Nsubs))]

        SubOff = np.concatenate(SubOff, axis=0)
        block.write(sstart,SubOff)
        #print('group %d done!'%(gend-1))

def init_offblock():
    dtype  = ('<i8', (6,))
    dsize  = dest_r['SubGroups/SubhaloLenType'].size
    nfile = dest_r['SubGroups/SubhaloLenType'].Nfile

    blockname = 'SubGroups/SubhaloOffsetType'
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
    block = init_offblock()
    comm.barrier()
    
    Ngroups = gend
    istart = Ngroups * rank // size
    iend = Ngroups * (rank+1) // size
    
    print('rank %d processing groups %d to %d'%(rank,istart,iend),flush=True)
    comm.barrier()
    


    write_subgroup_offset(block,istart,iend)
    
    