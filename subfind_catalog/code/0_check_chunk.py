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
        return False
    
    if not os.path.isfile(grpfile):
        print('*** %s does not exist!'%grpfile, flush=True)
        return False
        
    if "Subhalo" not in h5py.File(tabfile,'r').keys():
        print('*** %s does not have Subhalo!'%tabfile, flush=True)
        return False
    
    return True
        


        
    
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
        if cid in grp['PartType5']['ParticleIDs'][:]:
            found = 1
        return found
    

    found_minid = np.array([find_minid(i) for i in range(Nsubs)])
    trouble_list = (found_minid == 0).nonzero()[0]
    # pick out the ones we really care about
    sMass_trouble = sMass[trouble_list]
    mask_massive = ((sMass_trouble * 1e10 / 0.6774) > minsmass).nonzero()[0]
    massive_trouble_list = trouble_list[mask_massive]
    
    if len(trouble_list) > 0:
        print('Chunk %04d has %d subhalos, %d problematic, %d has stellar mass > %.1e'%(c, Nsubs, len(trouble_list), len(massive_trouble_list), minsmass), flush=True)
    
    
    if len(massive_trouble_list) > 0:
        #------------ Plot the largest problematic subhalo -------------------------
        isort = np.argsort(sMass[massive_trouble_list])[::-1]
        
        fig, ax = plt.subplots(1,2, figsize=(9,4.2))
        box = 250000.
        
        if len(massive_trouble_list) > 1:
            iplot = [0, 1]
        else:
            iplot = [0]
            
        for ii in iplot:
            itar = massive_trouble_list[isort[ii]]
            cid = cenID[itar]
            beg, end = sOff[itar], sOff[itar] + sLen[itar]

            minpos = tab['Subhalo']['SubhaloPos'][itar]
            group_idx = tab['Subhalo']['SubhaloGroupNr'][itar]
            pos1 = grp['PartType1']['Coordinates'][beg[1] : end[1]] - minpos
            pos4 = grp['PartType4']['Coordinates'][beg[4] : end[4]] - minpos

        

            pos1[pos1 > box/2] -= box
            pos1[pos1 < -box/2] += box
            pos4[pos4 > box/2] -= box
            pos4[pos4 < -box/2] += box

            mstar = sMass[itar]*1e10/0.6774

            ax[ii].scatter(pos1[:, 0], pos1[:, 1], label='DM', s=1)
            ax[ii].scatter(pos4[:, 0], pos4[:, 1], label='Star', s=1)

            ax[ii].scatter(0, 0, label='MinPos', color='cyan', s=50)
            ax[ii].set_title('C%d, G%d, S%d, Mstar %.1e'%(c, group_idx, itar, mstar))

        plt.savefig(savedir + '/trouble_chunk%d_sub%d_group%d.png'%(c, itar, group_idx), bbox_inches='tight')
        plt.close(fig)
    
    
    return trouble_list, sMass[trouble_list]*1e10/0.6774
        


def check_chunk(chunk_idxlist):
    # print('Processing chunk:',c,flush=True)
    cprob_list = []
    subid_list = []
    smass_list = []
    for c in chunk_idxlist:
        subdir  = subroot + '/chunk%d.%d/output/'%(c,maxgroup_list[c])
        tabfile = subdir + tabname
        grpfile = subdir + grpname

        data_complete = check_missing_data(tabfile, grpfile)

        if not data_complete:
            continue
        tab = h5py.File(tabfile, 'r')
        grp = h5py.File(grpfile, 'r')
        
        subid, smass = check_splitted_subhalo(c, tab, grp)
        subid_list.append(subid)
        cprob_list.append(np.ones_like(subid) * c)
        smass_list.append(smass)
        
    if len(subid_list) == 0:
        return None 
     
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
    parser.add_argument('--minsmass',default=0.,type=float,help='minimum stellar mass in Msun to report/plot, default 0')
    args = parser.parse_args()
    
    
    #--------------------------
    tabname = args.tabfile
    grpname = args.grpfile
    subroot = args.subroot
    cstart = int(args.cstart)
    savedir = args.savedir
    minsmass = args.minsmass
    
    if minsmass > 0:
        print("Will only report/plot problematic subhalos with stellar mass > %.1e Msun"%(minsmass), flush=True)

    #--------------create the savedir if needed-----------
    comm.barrier()
    if rank == 0:
        if not os.path.exists(savedir):
            os.mkdir(savedir)

    #-------------- Split tasks ------------------------------
    comm.barrier()
    chunk_list, maxgroup_list = get_subfind_chunk(subroot)

    if args.cend < cstart:
        cend = len(chunk_list)
    else:
        cend = int(args.cend)
    Nchunks = int(cend - cstart)
    
    chunk_num_perrank = Nchunks//size
    extra_chunk = Nchunks - (chunk_num_perrank * size)
    chunk_idxlist = np.arange(chunk_num_perrank) * size + rank + cstart
    print(f"chunk_num_perrank: {chunk_num_perrank}, extra_chunk {extra_chunk}.", flush=True)
    if(rank < extra_chunk):
        chunk_idxlist = np.hstack((chunk_idxlist, [chunk_num_perrank * size + rank + cstart]))
    
    istart = Nchunks * rank // size
    iend = Nchunks * (rank+1) // size
    
    if rank == 0:
        print('Total chunks:',Nchunks)

    comm.barrier()
    if(len(chunk_idxlist)>=2):
        print('Rank %03d will process %03d chunk %03d, %03d ... %03d'%(rank, len(chunk_idxlist), chunk_idxlist[0], chunk_idxlist[1], chunk_idxlist[-1]), flush=True)
    else:
        print('Rank %03d will process %03d chunk'%(rank, len(chunk_idxlist)), flush=True)
    comm.barrier()
    
    data = check_chunk(chunk_idxlist)    
    data_size = len(data)
    
    sendbuf_chunk = data['chunk'].flatten()
    sendbuf_subidx = data['subidx'].flatten()
    sendbuf_mstar = data['mstar'].flatten()
    
    comm.barrier()
    
    datasize_list = np.array(comm.gather(data_size, root=0))
    
    
    if rank == 0:
        print(f"rank {rank} get the total problematic data. Size: {datasize_list}")    
        recv_chunk = np.zeros(sum(datasize_list), dtype=np.int32)
        recv_subidx = np.zeros(sum(datasize_list), dtype=np.int32)
        recv_mstar = np.zeros(sum(datasize_list), dtype=np.float32)
    else:
        recv_chunk = None 
        recv_subidx = None 
        recv_mstar = None 

    comm.Gatherv(sendbuf_chunk, (recv_chunk, datasize_list), root=0)
    comm.Gatherv(sendbuf_subidx, (recv_subidx, datasize_list), root=0)
    comm.Gatherv(sendbuf_mstar, (recv_mstar, datasize_list), root=0)

    if rank == 0:
        print(f"rank {rank} get the gather data. Size: {len(recv_chunk)}, {len(recv_subidx)}, {len(recv_mstar)}.")
        
        dt = np.dtype([('chunk', np.int32), ('subidx', np.int32), ('mstar', np.float32)])
        gather_data = np.zeros(len(recv_chunk), dtype=dt)
        gather_data['chunk'] = recv_chunk
        gather_data['subidx'] = recv_subidx
        gather_data['mstar'] = recv_mstar
        if len(recv_chunk) > 0:
            np.save(os.path.join(savedir, "trouble_subhalo.npy"), gather_data)
    
            
            
            
            

        
        
        
        
        
        
        
    
    
    
    
    

        
        
        
                
            
            
                
        
        
        
        
        
        
        
        

        

















 
    
   
    
    


    
    
    
    
