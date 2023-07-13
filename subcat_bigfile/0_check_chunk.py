"""
A few completeness and consistency checks on the chunk subfind data
before we post-process them to get the bigfile subgroup catalog;
checks:
1. subfind is successfully run on all chunks (no missing data)
2. mpi issue of subgroup splitting (identified by minID not in subhalo)
3. count cases of large BH in small subgroup (todo)

feel free to add more
"""

import numpy as np
from bigfile import BigFile,FileMPI
import h5py
import sys,os
import glob
import argparse
from bf_util import *
from mpi4py import MPI
import matplotlib
import matplotlib.pyplot as plt
    

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


# in case the chunk data does not come with OffsetType
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
    


def check_missing_data(tabfile, grpfile):
    if not os.path.isfile(tabfile):
        print('*** %s does not exist!'%tabfile, flush=True)
    
    if not os.path.isfile(grpfile):
        print('*** %s does not exist!'%grpfile, flush=True)
        
    if "Subhalo" not in h5py.File(tabfile,'r').keys():
        print('*** %s does not have Subhalo!'%tabfile, flush=True)
        
        
    
def check_splitted_subhalo(c, tab, grp):
    sLen = tab['Subhalo']['SubhaloLenType'][:]
    sOff = tab['Subhalo']['SubhaloOffsetType'][:]
    cenID = tab['Subhalo']['SubhaloIDMostbound'][:]
    sMass = tab['Subhalo']['SubhaloMassType'][:][:,4]
    
    Nsubs = len(sLen)
    
    def find_minid(i):
        cid = cenID[i]
        beg, end = sOff[i], sOff[i] + sLen[i]
        found = 0
        for p in [0, 1, 4, 5]:
            if cid in grp['PartType%d'%p]['ParticleIDs'][beg[p] : end[p]]:
                found = 1
                break
        return found
    

    found_minid = np.array([find_minid(i) for i in range(Nsubs)])
    trouble_list = (found_minid == 0).nonzero()[0]
    
    
    if len(trouble_list) > 0:
        print('Chunk %04d has %d subhalos, %d problematic'%(c, Nsubs, len(trouble_list)), flush=True)
        
        #------------ Plot the largest problematic subhalo -------------------------
        ii = np.argmax(sMass[trouble_list])
        itar = trouble_list[ii] + 1
        cid = cenID[itar]
        beg, end = sOff[itar], sOff[itar] + sLen[itar]
        
        minpos = tab['Subhalo']['SubhaloPos'][itar]
        pos1 = grp['PartType1']['Coordinates'][beg[1] : end[1]] - minpos
        pos4 = grp['PartType4']['Coordinates'][beg[4] : end[4]] - minpos
        
        mstar = sMass[itar]*1e10/0.6774
        
        if mstar > 1e8:
            fig, ax = plt.subplots(1,1, figsize=(6,6))
            ax.scatter(pos1[:, 0], pos1[:, 1], label='DM', s=1)
            ax.scatter(pos4[:, 0], pos4[:, 1], label='Star', s=1)

            ax.scatter(0, 0, label='MinPos', color='cyan', s=50)
            ax.set_title('Chunk%d, Sub%d, Mstar %.1e'%(c, itar, mstar))

            plt.savefig(savedir + '/trouble_chunk%d_sub%d.png'%(c, itar), bbox_inches='tight')
            plt.close(fig)
    
    
    return trouble_list, sMass[trouble_list]*1e10/0.6774
        
            
def check_chunk(cstart, cend):
    # print('Processing chunk:',c,flush=True)
    cprob_list = []
    subid_list = []
    smass_list = []
    for c in range(cstart, cend):
        subdir  = subroot + '/chunk%d.%d/output/'%(c,maxgroup_list[c])
        tabfile = subdir + tabname
        grpfile = subdir + grpname

        check_missing_data(tabfile, grpfile)

        tab = h5py.File(tabfile, 'r')
        grp = h5py.File(grpfile, 'r')
        
        subid, smass = check_splitted_subhalo(c, tab, grp)
        subid_list.append(subid)
        cprob_list.append(np.ones_like(subid) * c)
        smass_list.append(smass)
        
    
     
    subid_list = np.concatenate(subid_list)
    cprob_list = np.concatenate(cprob_list)
    smass_list = np.concatenate(smass_list)
    
    dt = np.dtype([('chunk', np.int32), ('subidx', np.int32), ('mstar', np.float32)])
    data = np.zeros(len(subid_list), dtype=dt)
    
    data['chunk'] = cprob_list
    data['subidx'] = subid_list
    data['mstar'] = smass_list
    return data
        

#     for c in troubles.keys():
#         print('chunk %04d has %d problematic subhalos'%(c, len(troubles[c])))

        
        

if __name__ == "__main__":
    
    #----------- MPI Init ----------------
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    if rank == 0:
        print('Number of processes:', size, flush=True)

    #-------------- Cmd line Args ------------------------------
    parser = argparse.ArgumentParser(description='subgrpID-subfind')
    parser.add_argument('--subroot',required=True,type=str,help='path of the subfind directory')
    parser.add_argument('--tabfile',required=True,type=str,help='name of the tabfile')
    parser.add_argument('--grpfile',required=True,type=str,help='name of the subfind group file')
    parser.add_argument('--cstart',default=0,type=int,help='starting chunk')
    parser.add_argument('--cend',default=-1,type=int,help='ending chunk (exclusive)')
    parser.add_argument('--savedir',required=True,type=str,help='path for saving debug info')
    args = parser.parse_args()
    
    
    #--------------------------
    tabname = args.tabfile
    grpname = args.grpfile
    subroot = args.subroot
    cstart = int(args.cstart)
    savedir = args.savedir

    #-------------- Split tasks ------------------------------
    comm.barrier()
    chunk_list, maxgroup_list = get_subfind_chunk(subroot)

    if args.cend < cstart:
        cend = len(chunk_list)
    else:
        cend = int(args.cend)
    Nchunks = int(cend - cstart)
    
    istart = Nchunks * rank // size
    iend = Nchunks * (rank+1) // size
    
    if rank == 0:
        print('Total chunks:',Nchunks)

    comm.barrier()
    print('Rank %03d will process chunk %03d to chunk %03d'%(rank,cstart+istart,cstart+iend), flush=True)
    comm.barrier()
    
    chunks = chunk_list[cstart:cend]
    data = check_chunk(istart, iend)
    
            
            
            
            

        
        
        
        
        
        
        
    
    
    
    
    

        
        
        
                
            
            
                
        
        
        
        
        
        
        
        

        

















 
    
   
    
    


    
    
    
    
