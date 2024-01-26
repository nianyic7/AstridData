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
        sf_ID = []
        for g in grpfile:
            sbgrp = h5py.File(g,'r')
            # try except to deal with 0th chunk
            # where there is no 1 and 5
            try:
                sf_ID.append(sbgrp['PartType%d'%p]['ParticleIDs'][:])
            except KeyError:
                sf_ID.append(np.array([],dtype=np.int64))
        sf_ID = np.concatenate(sf_ID, axis=0)
        sub_sort[p] = np.argsort(sf_ID)
    return sub_sort



def subhalo_offset(tabfile):
    SubhaloLenType = np.concatenate([h5py.File(t,'r')['Subhalo']['SubhaloLenType'][:] for t in tabfile], axis=0)
    GroupLenType = np.concatenate([h5py.File(t,'r')['Group']['GroupLenType'][:] for t in tabfile], axis=0)
    GroupNsubs = np.concatenate([h5py.File(t,'r')['Group']['GroupNsubs'][:] for t in tabfile])
    GroupFirst = np.concatenate([h5py.File(t,'r')['Group']['GroupFirstSub'][:] for t in tabfile])
    
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
    
    # this may not exist for certain snapshots!
    SubhaloOffType = subhalo_offset(tabfile)
    SubhaloLenType = np.concatenate([h5py.File(t,'r')['Subhalo']['SubhaloLenType'][:] for t in tabfile], axis=0)
    SubhaloTotMass = np.concatenate([h5py.File(t,'r')['Subhalo']['SubhaloMass'][:] for t in tabfile], axis=0)

    # sort chunk subhalos by mass
    sort_mass = np.argsort(SubhaloTotMass)[::-1]
        
    grp_Nsubs      = np.concatenate([h5py.File(t,'r')['Group']['GroupNsubs'][:] for t in tabfile], axis=0)
    Npart          = np.sum([h5py.File(g,'r')['Header'].attrs['NumPart_ThisFile'] for g in grpfile], axis=0)
    
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


def reorder_group(ig,Length,Offset,sgpIDs):
    """
    Get the sorting index tp sort by subgroup IDs
    for particles in each 
    """

    # relative position in group
    start, end  = Offset[ig] - Offset[0], Offset[ig] - Offset[0] + Length[ig]
    # print('start,end',start,end)
    p = 1
    subID_set = set()
    for p in [0,1,4,5]:
        subID_set = subID_set.union(set(sgpIDs[p][start[p]:end[p]]))
    subID_arr   = np.array(list(subID_set))
    subID_arr   = np.sort(subID_arr)
    
    if len(subID_arr)==0 or subID_arr[0] == FUZZ:
        # no substructure in this group
        return None,None

    fof_sort = {}
    subID_sorted = {}
    for p in [0,1,4,5]:
        sub_ID      = sgpIDs[p][start[p]:end[p]]
        sort_ind    = np.argsort(sub_ID)
        fof_sort[p] = sort_ind
        subID_sorted[p] = sub_ID[sort_ind]
    return fof_sort,subID_sorted


def write_group_order(ig,Length,Offset,subID_sorted,new_order):
    start, end  = Offset[ig], Offset[ig] + Length[ig]
    
    for p in [0,1,4,5]:
    # for p in [1]:
        blockname1 = '%d/SubgroupID'%p
        blockname2 = '%d/NewIndex'%p
        if subID_sorted is None:
            data =  FUZZ * np.ones(end[p] - start[p])
            blocks[0][p].write(start[p], data)
            
            data = np.arange(end[p] - start[p])
            blocks[1][p].write(start[p], data)
        else:
            blocks[0][p].write(start[p], subID_sorted[p])
            blocks[1][p].write(start[p], new_order[p])


    
def process_chunk(c):
    print('Processing chunk:',c,flush=True)
    subdir  = subroot + '/chunk%d.%d/output/'%(c,maxgroup_list[c])
    if os.path.isdir(subdir + 'groups_%d'%(snap)):
        tabfile = sorted(glob.glob(subdir + 'groups_%03d/fof_subhalo_tab_%03d.*.hdf5'%(snap, snap)))
        grpfile = sorted(glob.glob(subdir + 'snapdir_%03d/snap_%03d.*.hdf5'%(snap, snap)))
    else:
        tabfile = [subdir + tab]
        grpfile = [subdir + grp]

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
    
    if "Subhalo" in h5py.File(tabfile[0],'r').keys():
        try: # assume arepo name
            subsChunk    = h5py.File(tabfile[0],'r')['Header'].attrs['Nsubgroups_Total']
        except KeyError: # revert to gadget4 naming
            subsChunk    = h5py.File(tabfile[0],'r')['Header'].attrs['Nsubhalos_Total']
        print('Total subgroups this chunk:',subsChunk, flush=True)
        #----------------- reorder and assign subID -----------------------
        fof_sort = fof_sorted(pstart, pend)
        sub_sort = sub_sorted(grpfile)
        sgpIDs, SubhaloLenType  = assign_sid_chunk(tabfile,grpfile)
        for p in [0,1,4,5]:
            sgpIDs[p][fof_sort[p]] = sgpIDs[p][sub_sort[p]]

        for ig in range(Ngroups):
            # i is the relative group index
            # gidx is the absolute group index
            new_order,subID_sorted = reorder_group(ig,Length,Offset,sgpIDs)
            write_group_order(ig,Length,Offset,subID_sorted,new_order)

    else: #no subhalo in the entire chunk
        #------------------ Nothing to do, just write -----------------------
        print('skipping entire chunk:',c,flush=True)
        for ig in range(Ngroups):
            # catagorize all as FUZZ
            # do not reorder groups in this chunk
            subID_sorted,new_order = None,None
            write_group_order(ig,Length,Offset,subID_sorted,new_order)

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
    parser.add_argument('--snap',required=True,type=int,help='snapshot number')
    parser.add_argument('--pigfile',required=True,type=str,help='path of the PIG file directory')
    parser.add_argument('--subroot',required=True,type=str,help='path of the subfind directory')
    parser.add_argument('--tabfile',required=True,type=str,help='name of the tabfile')
    parser.add_argument('--grpfile',required=True,type=str,help='name of the subfind group file')
    parser.add_argument('--dest',required=True,type=str,help='path of the output file directory')
    parser.add_argument('--cstart',default=0,type=int,help='starting chunk')
    
    args = parser.parse_args()
    
    #--------------------------
    pig = BigFile(args.pigfile)
    tab = args.tabfile
    grp = args.grpfile
    subroot = args.subroot
    cstart = int(args.cstart)

    snap = int(args.snap)
    
    dest  = FileMPI(comm, args.dest, create=True)
    gsize = pig['FOFGroups/LengthByType'].size
    FUZZ = 100000000
    


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

    
    Nchunks = int(len(chunk_list) - cstart)
    
    istart = Nchunks * rank // size
    iend = Nchunks * (rank+1) // size
    if rank == 0:
        print('Total chunks:',Nchunks)
        print('Saving dir:',args.dest,flush=True)
    comm.barrier()
    print('Rank %03d will process chunk %03d to chunk %03d'%(rank,istart,iend))
    comm.barrier()
    chunks = chunk_list[cstart:]
    for chunk in chunks[istart:iend]:
        process_chunk(chunk)
            
            
            
            

        
        
        
        
        
        
        
    
    
    
    
    

        
        
        
                
            
            
                
        
        
        
        
        
        
        
        

        

















 
    
   
    
    


    
    
    
    
