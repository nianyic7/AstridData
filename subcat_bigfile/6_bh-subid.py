import numpy as np
from bigfile import BigFile,FileMPI
import h5py
import sys,os
import glob
import argparse
from bf_util import *
from mpi4py import MPI
# from scipy.stats import rankdata

# """
# assign subgroupID to BHs
# also rewrite mbt cols using True BH mass

# """


#--------------------------------------------------------------------------
 
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
    # gend    = int(args.gend)

    
    comm.barrier()
    #---------- Initialize  blocks --------------
    blockname = '5-2/SubgroupIndex'
    dtype = 'i8'
    dsize = dest_r['5-2/ID'].size
    nfile = dest_r['5-2/ID'].Nfile
    if gstart == 0:
        block5 = dest_w.create(blockname, dtype, dsize, nfile)
    else:
        block5 = dest_w[blockname]
    comm.barrier()
    
#     blockname = '4/SubgroupIndex'
#     dtype = 'i8'
#     dsize = dest_r['4/ID'].size
#     nfile = dest_r['4/ID'].Nfile
#     if gstart == 0:
#         block4 = dest_w.create(blockname, dtype, dsize, nfile)
#     else:
#         block4 = dest_w[blockname]
#     comm.barrier()
        
        
    # ----------- Split tasks --------------------

    gLength  = dest_r['FOFGroups/LengthByType']
    gOffset  = dest_r['FOFGroups/OffsetByType']
    
    sLength  = dest_r['SubGroups/SubhaloLenType2']
    sOffset  = dest_r['SubGroups/SubhaloOffsetType2']
    
    FirstSub = dest_r['FOFGroups/GroupFirstSub'][:]
    Nsubs    = dest_r['FOFGroups/GroupNsubs']

    # number of subhalos with BHs
    NsubsTot   = sLength.size
    FinalSub   = (sLength[:][:,5] > 0).nonzero()[0][-1] + 1
    FinalGroup = (FirstSub >= FinalSub).nonzero()[0][0] + 1
    

    # starting and ending subgroup
    sstart = FinalSub * rank // size
    send   = FinalSub * (rank + 1) // size
    
    # starting and ending group
    istart = ((FirstSub >= 0) & (FirstSub <= sstart)).nonzero()[0][-1]
    iend   = (FirstSub >= send).nonzero()[0][0]
    
    

    if rank == 0:
        print('Sstart: %d, Total Subroups:'%sstart, FinalSub, "Total Groups:", FinalGroup)
        print('Saving dir:', args.dest,flush=True)
        lastBHinSub = sOffset[FinalSub][5]
        totalBHs = dest_r['5/ID'].size
        
        print('Total BHs:', totalBHs, 'Last BH in Subhalo:', lastBHinSub, flush=True)
        
    comm.barrier()
    print('Rank %d starts from Group %d to Group%d'%(rank,istart,iend),flush=True)
    comm.barrier()

    fsubs = FirstSub[istart : iend]
    nsubs = Nsubs[istart : iend]
    
    mask  = fsubs >= 0
    fsubs = fsubs[mask]
    nsubs = nsubs[mask]
    
    start, end = gOffset[istart][5], gOffset[iend][5]
    # minus one for default subgroup index
    sidx5 = - np.ones(end - start, dtype=np.int64)
    
    if len(fsubs) > 0:
        firstsub, lastsub = fsubs[0], fsubs[-1] + nsubs[-1]
        for s in range(firstsub, lastsub):
            sbeg, send = sOffset[s][5], sOffset[s][5] + sLength[s][5]
            if send <= sbeg:
                continue
            sidx5[sbeg : send] = s
    
    print('rank %d writing data from %d with length %d'%(rank, start, len(sidx5)), flush=True)
    block5.write(start, sidx5)
    print('rank %d done!'%rank,flush=True)
    
    comm.barrier()
    
    if rank == 0:
        print('writing the last BHs not in any subhalo...', flush=True)
        block5.write(lastBHinSub, - np.ones(totalBHs - lastBHinSub, dtype=np.int64))
        print('All done!', flush=True)
        
    
    

    #     #----------- Stars --------------------
    #     start, end = gOffset[istart][4], gOffset[iend][4]
    #     # minus one for default subgroup index
    #     sidx4 = - np.ones(end - start, dtype=np.int64)
    #     for s in range(firstsub, lastsub):
    #         sbeg, send = sOffset[s][4], sOffset[s][4] + sLength[s][4]
    #         if sbeg <= send:
    #             continue
    #         sidx4[sbeg : send] = s

    #     block4.write(start, sidx4)


        



