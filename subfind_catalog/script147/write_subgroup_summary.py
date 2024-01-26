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


def get_GroupSidx(gidx):
    start,end = Offset[gidx], Offset[gidx] + Length[gidx]
    # DO NOT assume that each subgroup has a DM particle!!
    subID = np.concatenate([dest_r['%d/SubgroupID'%p][start[p]:end[p]] for p in [0,1,4,5]])

    subID_set   = set(subID)
    if len(subID_set) == 0 or min(list(subID_set)) == FUZZ:
        # no substructure in this group
        return np.array([],dtype=np.int32)
    
    if max(list(subID_set)) == FUZZ:
        # there is inner fuzz
        subID_set.remove(FUZZ)
    subID_arr = np.sort(np.array(list(subID_set)))
    return subID_arr



def get_SubgroupOff(GroupOff,SubLen):
    if len(SubLen) == 0:
        return np.zeros((0,6),dtype = np.int32)
    SubOff = np.zeros_like(SubLen,dtype=int)
    SubOff[1:] = np.cumsum(SubLen[:-1],axis=0)
    return SubOff + GroupOff


def write_subgroups_chunk(ChunkFirstSub,SidxFOF2Chunk,GroupSidxChunk,SubOff,Subhalo):
    if len(GroupSidxChunk) == 0:
        return 
    assert len(GroupSidxChunk) == len(SubOff), \
    'Calculated offset does not match subgroup number! %d %d'%(len(GroupSidxChunk), len(SubOff))
    
    assert len(GroupSidxChunk) == len(SubOff), \
    'Calculated offset does not match subgroup number! %d %d'%(len(GroupSidxChunk), len(SubOff))
    
    new_ind = SidxFOF2Chunk[GroupSidxChunk]
    for block,name in blocklists:
        if name != 'SubGroups/SubhaloOffsetType':
            # reorder group particles
            feature = name.split('/')[1]
            data    = np.concatenate([S[feature][:] for S in Subhalo], axis=0)[new_ind]
            # write to file
            block.write(ChunkFirstSub,data)
        else:
            # need to recompute this field
            print('writing offset type starting with',SubOff[0],flush=True)
            block.write(ChunkFirstSub,SubOff)
        

def write_groups_chunk(gstart,FirstSubList,NsubsList):
    print('writing group:',gstart)
    blocklists_fof[0].write(gstart,np.array(FirstSubList))
    blocklists_fof[1].write(gstart,np.array(NsubsList))
    
    

def process_chunk(c, ChunkFirstSub):
    """
    write summary statistics for each chunk
    write Nsubs, first sub for each group
    """
    # groups in this chunk
    print('Processing chunk:',c,flush=True)
    subdir  = subroot + '/chunk%d.%d/output/'%(c,maxgroup_list[c])
    
    if os.path.isdir(subdir + 'groups_%d'%(snap)):
        tabfile = sorted(glob.glob(subdir + 'groups_%03d/fof_subhalo_tab_%03d.*.hdf5'%(snap, snap)))
    else:
        tabfile = [subdir + tab]

    if c == 0:
        gstart, gend = 0, maxgroup_list[c]
    else:
        gstart, gend = maxgroup_list[c-1], maxgroup_list[c]
        
    if "Subhalo" in h5py.File(tabfile[0],'r').keys():
        Subhalo = [h5py.File(t,'r')['Subhalo'] for t in tabfile]
        SubhaloTotMass = np.concatenate([S['SubhaloMass'][:] for S in Subhalo], axis=0)
        SubLengthChunk = np.concatenate([S['SubhaloLenType'][:] for S in Subhalo], axis=0)
    else:
        # no subhalo in the entire chunk
        # nothing to do
        print('No subgroups at all in chunk%d'%c,flush=True)
        return  ChunkFirstSub
        
    SidxFOF2Chunk  = np.argsort(SubhaloTotMass)[::-1]
    GroupSidxList  = [get_GroupSidx(gidx) for gidx in range(gstart,gend)]
    GroupSidxChunk = np.concatenate(GroupSidxList)
    
    SubLenList     = [SubLengthChunk[SidxFOF2Chunk[s]] for s in GroupSidxList]
    NsubsChunk     = np.array([len(s) for s in GroupSidxList])
    NsubsTot       = sum(NsubsChunk)
    
    FirstSubs      = np.zeros_like(NsubsChunk) + ChunkFirstSub
    FirstSubs[1:]  +=  np.cumsum(NsubsChunk[:-1])
    FirstSubs[NsubsChunk == 0] = -1
    
    GroupOff       = Offset[gstart:gend]
    
    # print('groupoff for group %d to %d'%(gstart,gend),GroupOff,flush=True)
    SubOffChunk    = [get_SubgroupOff(GroupOff[i],SubLenList[i]) for i in range(len(GroupOff))]
    SubOffChunk    = np.concatenate(SubOffChunk,axis=0).astype(int)
#     assert len(SubLengthChunk)==NsubsTot, \
#     "Nsub chunk does not match: fof=%d, chunk=%d"%(len(SubLengthChunk),NsubsTot)
    
#     print('sgroupoff for group %d to %d'%(gstart,gend),SubOffChunk,flush=True)
    
    write_groups_chunk(gstart,FirstSubs,NsubsChunk)
    write_subgroups_chunk(ChunkFirstSub,SidxFOF2Chunk,GroupSidxChunk,SubOffChunk,Subhalo)


    
def get_Nsubs(chunks,istart,iend):
    Nsubs = []
    for c in chunks[istart:iend]:
        subdir  = subroot + '/chunk%d.%d/output/'%(c,maxgroup_list[c])
        if os.path.isdir(subdir + 'groups_%d'%(snap)):
            tabfile = sorted(glob.glob(subdir + 'groups_%03d/fof_subhalo_tab_%03d.*.hdf5'%(snap, snap)))
        else:
            tabfile = [subdir + tab]
        if "Subhalo" in h5py.File(tabfile[0],'r').keys():
            try:
                NsubsChunk  = h5py.File(tabfile[0],'r')['Header'].attrs['Nsubgroups_Total']
            except KeyError:
                NsubsChunk  = h5py.File(tabfile[0],'r')['Header'].attrs['Nsubhalos_Total']
        else:
            NsubsChunk  = 0
        Nsubs.append(NsubsChunk)
    return Nsubs



def init_subblocks(NSubs):
    tabfile = glob.glob(subroot + '/chunk1.*/output/' + tab)[0]
    Subhalo = h5py.File(tabfile,'r')['Subhalo']
    
    blocklists = []
    for k in Subhalo.keys():
        blockname = 'SubGroups/%s'%k
        
        if cstart == 0:
            dtype  = Subhalo[k].dtype
            dshape = Subhalo[k].shape
            dsize = NSubs
            nfile = 256 # same as FOFGroups
            if len(dshape) > 1:
                width = dshape[1]
                dtype = np.dtype((dtype,(width,)))

            block = dest_w.create(blockname, dtype, dsize, nfile)
            if rank == 0:
                print('created block %s with dsize %d'%(blockname,dsize),flush=True)
        else:
            block = dest_w[blockname]
            if rank == 0:
                print('block %s exists, dsize=%d'%(blockname,block.size),flush=True)
            
        blocklists.append((block,blockname))

    if 'SubhaloOffsetType' not in Subhalo.keys():
        # init offset block if there isn't one
        blockname = 'SubGroups/SubhaloOffsetType'
        if cstart == 0:
            dtype  = ('i8',(6,))
            dshape = Subhalo['SubGroups/SubhaloLenType'].shape
            dsize = NSubs
            nfile = 256 # same as FOFGroups
            block = dest_w.create(blockname, dtype, dsize, nfile)
        else:
            block = dest_w[blockname]
        blocklists.append((block,blockname))

    # New FOFGroup Blocks
    dtype = np.dtype('i4')
    dsize = pig['FOFGroups/GroupID'].size
    nfile = pig['FOFGroups/GroupID'].Nfile
    blocklists_fof = []
    for feature in ['GroupFirstSub','GroupNsubs']:
        blockname = 'FOFGroups/%s'%feature
        if cstart == 0:
            block = dest_w.create(blockname, dtype, dsize, nfile)
        else:
            block = dest_w[blockname]
        blocklists_fof.append(block)
        if rank == 0:
            print('Initizlized block:', blockname,flush=True)
    
    return blocklists,blocklists_fof

#--------------------------------------------------------------------------
 
if __name__ == "__main__":
    #----------- MPI Init ----------------
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    #-------------- Cmd line Args ------------------------------
    parser = argparse.ArgumentParser(description='subgrpID-subfind')
    parser.add_argument('--pigfile',required=True,type=str,help='path of the PIG file directory')
    parser.add_argument('--snap',required=True,type=int,help='snapshot number')
    parser.add_argument('--subroot',required=True,type=str,help='path of the subfind directory')
    parser.add_argument('--tabfile',required=True,type=str,help='name of the tabfile')
    parser.add_argument('--grpfile',required=True,type=str,help='name of the subfind group file')
    parser.add_argument('--dest',required=True,type=str,help='path of the output file directory')
    parser.add_argument('--cstart',default=0,type=int,help='starting chunk')
    parser.add_argument('--minpart',default=200,type=int,help='min particle in halo to start rewriting')

    args = parser.parse_args()
    
    #--------------------------
    pig = BigFile(args.pigfile)
    tab = args.tabfile
    grp = args.grpfile
    subroot = args.subroot
    cstart = int(args.cstart)
    minpart = int(args.minpart)
    snap = int(args.snap)
    
    dest_w  = FileMPI(comm, args.dest, create=True)
    dest_r  = BigFile(args.dest)
    gsize   = pig['FOFGroups/LengthByType'].size
    FUZZ    = 100000000
    
    Length  = pig['FOFGroups/LengthByType']
    Offset  = pig['FOFGroups/OffsetByType']

    gend = 0
    if rank == 0:
        LengthByType = Length[:]
        N01 = (LengthByType[:,0] + LengthByType[:,1]).astype(int)
        gend = np.searchsorted(-N01, -minpart)
        del LengthByType, N01
        print('Ending Gidx for reordering: %d'%gend,flush=True)
        
    comm.barrier()
        
    gend = comm.allreduce(gend)
    if rank == 1:
        print('Ending Gidx for reordering: %d'%gend,flush=True)
    comm.barrier()
    
    # ----------- Split tasks --------------------
    chunk_list, maxgroup_list = get_subfind_chunk(subroot)
    cend = group_chunk_dir(gend) + 1
    # cend = len(chunk_list)
    
    Nchunks = int(cend - cstart)
    
    Nproc1 = max(1,size//5)
    Nproc2 = size - Nproc1
    
    batch2 = Nchunks//2
    batch1 = Nchunks - batch2
    
    if rank <  Nproc1:
        istart = batch1 * rank // Nproc1
        iend = batch1 * (rank+1) // Nproc1
    else:
        istart = batch1 + batch2 * (rank-Nproc1) // Nproc2
        iend = batch1 + batch2 * (rank-Nproc1+1) // Nproc2
        
    chunks = chunk_list[cstart:cend]
    
    #-------------------- Check -----------------------
    nchunks = iend - istart
    nchunks_all = comm.allreduce(nchunks)
    assert nchunks_all==Nchunks, \
    'tasks not splitted correctly, Should have %d, splitted %d'%(Nchunks,nchunks_all)
    #--------------------------------------------------
        
        
    if rank == 0:
        print('Total chunks:',Nchunks,flush=True)
        print('Saving dir:',args.dest,flush=True)
    comm.barrier()
    print('Rank %03d will process chunk %03d to chunk %03d'%(rank,istart,iend),flush=True)
    comm.barrier()
    

    #------------ First pass to figure out total data size -----------------
    NsubsList = get_Nsubs(chunks,istart,iend)
    print('Nsubslist for rank %d'%rank,istart,iend,NsubsList,flush=True)
    Nsubs = sum(NsubsList)
    # how many subgroups there are before this task
    Nsubs_pre = sum(get_Nsubs(chunk_list,0,cstart))
    
    comm.barrier()
    TaskOff   = sum(comm.allgather(Nsubs)[:rank]) + Nsubs_pre
    # total subgroups
    NsubsTot = comm.allreduce(Nsubs)
    comm.barrier()
    if rank == 0:
        print(NsubsTot,flush=True)
        print('Nsubs before cstart: %d'%Nsubs_pre)
        print('Max chunk %d, maxgroup %d, Total Subgroups: %d'%(cend,maxgroup_list[cend-1],NsubsTot),flush=True)
    comm.barrier()
    
    print('Offset for rank %03d: %d'%(rank,TaskOff))
    print('Subgroup Length for rank %03d: %d'%(rank,Nsubs))
    comm.barrier()
    
    
    #---------- Initialize all blocks --------------
    blocklists,blocklists_fof = init_subblocks(NsubsTot)
    comm.barrier()
    
    #---------- Do the Work --------------
    FirstSub = TaskOff
    for i,chunk in enumerate(chunks[istart:iend]):
        NsubsChunk = NsubsList[i]
        process_chunk(chunk,FirstSub)
        FirstSub += NsubsChunk
        print('chunk %d done! firstsub=%d'%(chunk,FirstSub),flush=True)
    
    