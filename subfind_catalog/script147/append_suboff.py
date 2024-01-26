import numpy as np
from bigfile import BigFile,FileMPI
import h5py
import sys,os
import glob
import argparse
from bf_util import *
from mpi4py import MPI



def get_SubgroupOff(GroupOff,SubLen):
    if len(SubLen) == 0:
        return np.zeros((0,6),dtype = np.int64)
    SubOff = np.zeros_like(SubLen,dtype=np.int64)
    SubOff[1:] = np.cumsum(SubLen[:-1],axis=0)
    return (SubOff + GroupOff).astype(np.int64)


def write_subgroup_offset(block,gstart,gend):
    if gstart < 1000:
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
        # print('group %d done!'%(gend-1))

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
    parser.add_argument('--gend',required=True,type=int,help='min particle in halo to start rewriting')

    args = parser.parse_args()
    
    #--------------------------
    gend = int(args.gend)
    dest_w  = FileMPI(comm, args.dest, create=True)
    dest_r  = BigFile(args.dest)
    gsize   = dest_r['FOFGroups/LengthByType'].size
    ssize   = dest_r['SubGroups/SubhaloLenType'].size
    FUZZ    = 100000000
    
    Length  = dest_r['FOFGroups/LengthByType']
    Offset  = dest_r['FOFGroups/OffsetByType']

    #---------- Initialize all blocks --------------
    block = init_offblock()
    comm.barrier()
    
    Ngroups = gend
    istart = Ngroups * rank // size
    iend = Ngroups * (rank+1) // size
    
    print('rank %d processing groups %d to %d'%(rank,istart,iend),flush=True)
    comm.barrier()
    


    write_subgroup_offset(block,istart,iend)
    
    