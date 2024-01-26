import numpy as np
from bigfile import BigFile,FileMPI
import h5py
import sys,os
import glob
import argparse
from bf_util import *
from mpi4py import MPI



def rewrite_group(igroup,blockname,p):
    start,end = Starts[igroup],  Ends[igroup]
    NewOrder  = dest_r['%d/NewIndex'%p][start:end]
    data      = pig[blockname][start:end][NewOrder]
    block.write(start,data)



def init_blocks(source,dest,blocknames):
    pig = BigFile(source)
    blocklists = {p:[] for p in [0,1,4,5]}
    for blockname in blocknames:
        p = int(blockname.split('/')[0])
        dtype,dsize,nfile = dtype_size_nfile(pig,blockname)
        if gstart == 0:
            block = dest.create(blockname, dtype, dsize, nfile)
        else:
            block = dest[blockname]
        if rank == 0:
            print('Initizlized block:', blockname,flush=True)
        blocklists[p].append((block, blockname))
    return blocklists

#--------------------------------------------------------------------------
 
if __name__ == "__main__":
    #----------- MPI Init ----------------
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()


    #-------------- Cmd line Args ------------------------------
    parser = argparse.ArgumentParser(description='subgrpID-subfind')
    parser.add_argument('--pigfile',required=True,type=str,help='path of the PIG file directory')
    parser.add_argument('--dest',required=True,type=str,help='path of the output file directory')
    parser.add_argument('--gstart',default=0,type=int,help='where to begin the rewriting')
    parser.add_argument('--gend',default=0,type=int,help='where to finish the rewriting')
    parser.add_argument('--minpart',default=200,type=int,help='min particle in halo to start rewriting')
    parser.add_argument('--blocknames', nargs='+', help='blocks to process', required=True)
    
    args = parser.parse_args()
    
    #--------------------------
    pig = BigFile(args.pigfile)
    
    dest_w  = FileMPI(comm, args.dest, create=True)
    dest_r  = BigFile(args.dest)
    gsize   = pig['FOFGroups/LengthByType'].size
    gstart  = int(args.gstart)
    gend = int(args.gend)
    
    minpart = int(args.minpart)
    blocknames = list(args.blocknames)
    if rank == 0:
        print('blocks to process:', blocknames,flush=True)
    
    
    comm.barrier()
    #---------- Initialize all blocks --------------
    blocklists = init_blocks(source=args.pigfile, dest=dest_w, blocknames=blocknames)
    comm.barrier()
        
    Length  = pig['FOFGroups/LengthByType']
    Offset  = pig['FOFGroups/OffsetByType']
    
    gfinal = 0
    if rank == 0:
        LengthByType = Length[:]
        N01 = (LengthByType[:,0] + LengthByType[:,1]).astype(int)
        gfinal = np.searchsorted(-N01, -minpart)
        del LengthByType, N01
        print('Ending Gidx for reordering: %d'%gfinal,flush=True)
        
    comm.barrier()
        
    gfinal = comm.allreduce(gfinal)
    if rank == 1:
        print('Abs Final Gidx for reordering: %d'%gfinal,flush=True)
    comm.barrier()
    
    # cap at minpart thresh group

    if gfinal < int(gend):
        gend = gfinal
        
    
    # ----------- Split tasks --------------------
    Ngroups = gend - gstart
    
    Nproc2 = max(1,size//2)
    Nproc1 = size - Nproc2
    print("Proc1: %d Proc2: %d"%(Nproc1,Nproc2),flush=True)
    
    batch1 = max(1000, Ngroups//100)
    batch2 = Ngroups - batch1

    print('Batch1: %d Batch2: %d'%(batch1,batch2),flush=True)
    
    if rank <  Nproc1:
        istart = batch1 * rank // Nproc1 + gstart
        iend = batch1 * (rank + 1) // Nproc1 + gstart
    else:
        istart = batch1 + batch2 * (rank-Nproc1) // Nproc2 + gstart
        iend = batch1 + batch2 * (rank-Nproc1+1) // Nproc2 + gstart

    print('Rank %d starts from Group %d to Group%d'%(rank,istart,iend),flush=True)
    if rank == 0:
        print('Gstart: %d Gend: %d Total groups:'%(gstart,gend),Ngroups)
        print('Saving dir:',args.dest,flush=True)
    comm.barrier()

    
    #--------------------------------------------------------------


    for p in [4,5,0,1]:
        if len(blocklists[p]) == 0:
            continue
        #------------------------------------------------
        Starts = Offset[istart:iend][:,p]
        Lens   = Length[istart:iend][:,p]
        Ends   = Starts + Lens

        if (gstart <= 1000000) and (istart <= max(Ngroups//20,300000)):
            # work on groups one by one
            for i,gidx in enumerate(np.arange(istart,iend)):
                for block,name in blocklists[p]:
                    rewrite_group(i,name,p)
            print('Rank %d finished all'%(rank),flush=True)
        else:
            # bundle data
            off       = Starts - Starts[0:1]
            offset    = np.concatenate([np.ones(Lens[i]) * off[i] for i in range(iend-istart)]).astype(int)
            del off
            NewOrder  = dest_r['%d/NewIndex'%p][Starts[0]:Ends[-1]] + offset
            del offset
            for block,name in blocklists[p]:
                data      = pig[name][Starts[0]:Ends[-1]][NewOrder]
                block.write(Starts[0],data)
            print('Rank %d finished all'%(rank),flush=True)
            
            
        #---------------- copy over the smaller groups ----------------
        if (gend==gfinal) and (gfinal < gsize):
            if rank == size//2:
                rstart = Offset[gfinal]
                print('rank %d copying over the small groups starting at:'%rank, rstart)
                for block,name in blocklists[p]:
                    data = pig[name][rstart[p]:]
                    block.write(rstart[p],data)
                



