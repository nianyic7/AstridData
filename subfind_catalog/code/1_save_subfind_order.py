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


def fof_sorted(start,end):
    fof_sort = {}
    for p in [0,1,4,5]:
        fof_ID      = pig['%d/ID'%p][start[p]:end[p]]
        fof_idx     = np.argsort(fof_ID)
        fof_sort[p] = fof_idx
    return fof_sort


def sub_sorted(grpfile):
    sub_sort = {}
    for p in [0,1,4,5]:
        # load sf_ID* and Group, Subhalo offset from subfind
        sbgrp = h5py.File(grpfile,'r')
        # try except to deal with 0th chunk
        # where there is no 1 and 5
        try:
            sf_ID = sbgrp['PartType%d'%p]['ParticleIDs'][:]
        except KeyError:
            sf_ID = np.array([],dtype=np.int64)
        sub_sort[p] = np.argsort(sf_ID)
    return sub_sort



def subhalo_offset(tabfile):
    SubhaloLenType = h5py.File(tabfile,'r')['Subhalo']['SubhaloLenType'][:]
    GroupLenType = h5py.File(tabfile,'r')['Group']['GroupLenType'][:]
    GroupNsubs = h5py.File(tabfile,'r')['Group']['GroupNsubs'][:]
    GroupFirst = h5py.File(tabfile,'r')['Group']['GroupFirstSub'][:]
    
    SubhaloOffsetType = np.zeros_like(SubhaloLenType,dtype=int)
    GroupOff = np.zeros(6,dtype=int)
    SubOff = np.zeros(6,dtype=int)
    for ig,lgroup in enumerate(GroupLenType):
        SubOff[:] = GroupOff[:]
        nsub = GroupNsubs[ig]
        firstsub = GroupFirst[ig]
        if firstsub == -1:
            GroupOff += lgroup
            continue
        for isub in range(firstsub,firstsub + nsub):
            SubhaloOffsetType[isub] = SubOff
            SubOff += SubhaloLenType[isub]
        GroupOff += lgroup
    return SubhaloOffsetType
    



def assign_sid_chunk(tabfile,grpfile):
    """
    returns: 
    
    sgpIDs : dict with sgpIDs[p] the subgroupID of all particles type p, in chunk ordering
    SubhaloLenType : subhalo lbt of chunk data
    
    """
    
    # this may not exist for certain snapshots!
    SubhaloOffType = subhalo_offset(tabfile)
    SubhaloLenType = h5py.File(tabfile,'r')['Subhalo']['SubhaloLenType'][:]
    SubhaloTotMass = h5py.File(tabfile,'r')['Subhalo']['SubhaloMass'][:]

    # sort chunk subhalos by mass
    sort_mass = np.argsort(SubhaloTotMass)[::-1]
        
    grp_Nsubs      = h5py.File(tabfile,'r')['Group']['GroupNsubs'][:]
    Npart          = h5py.File(grpfile,'r')['Header'].attrs['NumPart_ThisFile']
    
    Ngroups = len(grp_Nsubs)
    Nsubs   = len(SubhaloLenType)
    # --------------------------------------------------------------------
    # if we only care about subgroups then
    # we can assign unique subgroup ids in this chunk?
    sgpIDs = {}
    for p in [0,1,4,5]:
        # default to large integer so that the wandering 
        # particles are sorted to the end
        # better re-order all subgroups in this chunk by mass
        # to guarantee that old fof subgroups following the 
        # mass ordering
        sf_sgpID = FUZZ * np.ones(Npart[p],dtype='i4')
        scount = 0
        for i,s in enumerate(sort_mass):
            sN = SubhaloLenType[s]
            sO = SubhaloOffType[s]
            sf_sgpID[sO[p]:sO[p]+sN[p]]=i
        sgpIDs[p] = sf_sgpID
    return sgpIDs, SubhaloLenType



def get_group_newidx(p, ig, start, end, sgpIDs):
    """
    Get the sorting index and sorted SubIDs 
    (needed in case a group is splitted)
    for particles in each group
    
    Input:
    ig : relative group idx in this batch
    p: particle type
    Length: lbt of groups in this batch
    Offset: obt of groups in this batch
    sgpIds: dict of unsorted sID
    """

    # subID_arr = np.concatenate([sgpIDs[p][start[p]:end[p]] for p in [0,1,4,5]])
    if start[p] == end[p]:
        # for some reason no data...
        sort_ind = np.array([], dtype=int)
    elif np.min(sgpIDs[p][start[p]:end[p]]) == FUZZ:
        # no substructure in this group
        sort_ind = np.arange(end[p] - start[p]).astype(int)
    else:
        sort_ind    = np.argsort(sgpIDs[p][start[p]:end[p]]).astype(int)
    return sort_ind
    
    
    
    
def get_chunk_newidx_sid(p, Length, Offset, sgpIDs, Ngroups, skipchunk=False):
    starts, ends = (Offset - Offset[0:1]).astype(int), (Offset - Offset[0:1] + Length).astype(int)
    print("starting shape:", starts.shape, ends.shape, flush=True)
    
    if skipchunk:
        chunk_newidx = [np.arange(ends[ig][p] - starts[ig][p]) for ig in range(Ngroups)]
        chunk_sids = FUZZ * np.ones(ends[-1][p] - starts[0][p])
    else:
        chunk_newidx = [get_group_newidx(p, ig, starts[ig], ends[ig], sgpIDs) for ig in range(Ngroups)]

        chunk_sids = np.concatenate([sgpIDs[p][starts[ig][p]:ends[ig][p]][sort_ind] \
                                     for ig,sort_ind in enumerate(chunk_newidx) if starts[ig][p] < ends[ig][p]])


    return np.concatenate(chunk_newidx), chunk_sids
    
    
def write_chunk_data(p, start, end, subID_sorted,new_order):
    blocks[0][p].write(start[p], subID_sorted)
    blocks[1][p].write(start[p], new_order)

    
def process_chunk(c):
    print('Processing chunk:',c,flush=True)
    subdir  = subroot + '/chunk%d.%d/output/'%(c,maxgroup_list[c])
    tabfile = subdir + tab
    grpfile = subdir + grp

    if c == 0:
        gstart, gend = 0, maxgroup_list[c]
    else:
        gstart, gend = maxgroup_list[c-1], maxgroup_list[c]
    print('gstart',gstart,'gend',gend, flush=True)
        
    Length  = pig['FOFGroups/LengthByType'][gstart:gend]
    Offset  = pig['FOFGroups/OffsetByType'][gstart:gend]
    Ngroups = len(Length)
    
    # starting and ending particle index of this chunk in FOFGroups
    pstart, pend = Offset[0], Offset[-1] + Length[-1]
    
    if "Subhalo" in h5py.File(tabfile,'r').keys():
        try: # assume arepo name
            subsChunk    = h5py.File(tabfile,'r')['Header'].attrs['Nsubgroups_Total']
        except KeyError: # revert to gadget4 naming
            subsChunk    = h5py.File(tabfile,'r')['Header'].attrs['Nsubhalos_Total']
        print('Total subgroups this chunk:',subsChunk, flush=True)
        #----------------- reorder and assign subID -----------------------
        fof_sort = fof_sorted(pstart, pend)
        sub_sort = sub_sorted(grpfile)
        
        sgpIDs, SubhaloLenType  = assign_sid_chunk(tabfile,grpfile)
        
        for p in [0,1,4,5]:
            sgpIDs[p][fof_sort[p]] = sgpIDs[p][sub_sort[p]]
            
        for p in [0,1,4,5]:
            # no this type of particle in this chunk, nothing to write. 
            if pstart[p] == pend[p]:
                continue
            # get orders
            new_order, subID_sorted = get_chunk_newidx_sid(p, Length, Offset, sgpIDs, Ngroups)
            # write orders
            write_chunk_data(p, pstart, pend, subID_sorted, new_order)


    else: #no subhalo in the entire chunk
        #------------------ Nothing to do, just write -----------------------
        print('skipping entire chunk:',c,flush=True)
        sgpIDs = None
        for p in [0,1,4,5]:
            # no this type of particle in this chunk, nothing to write. 
            if pstart[p] == pend[p]:
                continue
            # get orders
            new_order, subID_sorted = get_chunk_newidx_sid(p, Length, Offset, sgpIDs, Ngroups, skipchunk=True)
            # write orders
            write_chunk_data(p, pstart, pend, subID_sorted, new_order)

    print('chunk %d done!'%c,flush=True)



if __name__ == "__main__":
    
    #----------- MPI Init ----------------
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    print('Number of processes:', size, flush=True)
    print('This task:', rank, flush=True)


    #-------------- Cmd line Args ------------------------------
    parser = argparse.ArgumentParser(description='subgrpID-subfind')
    parser.add_argument('--pigfile',required=True,type=str,help='path of the PIG file directory')
    parser.add_argument('--subroot',required=True,type=str,help='path of the subfind directory')
    parser.add_argument('--tabfile',required=True,type=str,help='name of the tabfile')
    parser.add_argument('--grpfile',required=True,type=str,help='name of the subfind group file')
    parser.add_argument('--dest',required=True,type=str,help='path of the output file directory')
    parser.add_argument('--cstart',default=0,type=int,help='starting chunk')
    parser.add_argument('--cend',default=-1,type=int,help='ending chunk (exclusive)')
    args = parser.parse_args()
    
    #--------------------------
    pig = BigFile(args.pigfile)
    tab = args.tabfile
    grp = args.grpfile
    subroot = args.subroot
    cstart = int(args.cstart)

    
    dest  = FileMPI(comm, args.dest, create=True)
    gsize = pig['FOFGroups/LengthByType'].size
    FUZZ = 100000000
    
    if rank==0:
        print('subfind root:', subroot, flush=True)

#    ------------------ Initialize Blocks --------------------
    comm.barrier()
    blocks = [{},{}]
    for p in [0,1,4,5]:
        dtype,dsize,nfile = dtype_size_nfile(pig,blockname='%d/ID'%(p))

        dtype = 'i4'
        blockname = '%d/SubgroupID'%(p)
        if cstart == 0:
            block = dest.create(blockname,dtype, dsize, nfile)
        else:
            block = dest[blockname]
        blocks[0][p] = block
        if rank == 0:
            print('Initialized block:', blockname,flush=True)

        blockname = '%d/NewIndex'%(p)
        if cstart == 0:
            block = dest.create(blockname,dtype, dsize, nfile)
        else:
            block = dest[blockname]
        blocks[1][p] = block
        if rank == 0:
            print('Initialized block:', blockname,flush=True)


    #-------------- Split tasks ------------------------------
    comm.barrier()
    chunk_list, maxgroup_list = get_subfind_chunk(subroot)

    if args.cend <= cstart:
        cend = len(chunk_list)
    else:
        cend = int(args.cend)
    Nchunks = int(cend - cstart)
    
    istart = Nchunks * rank // size
    iend = Nchunks * (rank+1) // size
    
    if rank == 0:
        print('Total chunks:',Nchunks)
        print('Saving dir:',args.dest,flush=True)
    comm.barrier()
    print('Rank %03d will process chunk %03d to chunk %03d'%(rank,cstart+istart,cstart+iend))
    comm.barrier()
    chunks = chunk_list[cstart:cend]
    for chunk in chunks[istart:iend]:
        process_chunk(chunk)
            
            
            
            

        
        
        
        
        
        
        
    
    
    
    
    

        
        
        
                
            
            
                
        
        
        
        
        
        
        
        

        

















 
    
   
    
    


    
    
    
    
