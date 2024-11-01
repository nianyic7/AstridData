import numpy as np
from bigfile import BigFile
import glob,os,struct
import pickle
import warnings
import argparse
import time

hh=0.6774

#-------------------------------------------------------

def load_bh_sft(sprev,sbeg):
    pig = BigFile('/hildafs/datasets/Asterix/PIG2/PIG_%03d'%sprev)
    battr = pig["Header"].attrs
    amin = battr["Time"][0]
    print('Minimum seeding atime:',amin,flush=True)
    
    bh_sft = {}
    for ff in sorted(glob.glob('/hildafs/datasets/Asterix/PIG2/PIG_*')):
        if ff[-3:]=='ind':
            continue
        if int(ff[-3:]) < sbeg:
            continue
        pig = BigFile(ff)
        print('Loading:',ff[-3:],flush=True)

        form = pig.open('5/StarFormationTime')[:]
        mask = form >= amin
        print('Num BHs in this file:',len(form[mask]),flush=True)
        bhid = pig.open('5/ID')[:][mask]
        form = form[mask]

        bh_sft.update({b:form[i] for i,b in enumerate(bhid)})
        print('Total Num BHs:',len(bh_sft),flush=True)
    return bh_sft
    

def assign_chunks(cbeg,csize,bh_sft):
    istart = cbeg
    bhid = np.fromiter(bh_sft.keys(), dtype=np.int64)
    form = np.fromiter(bh_sft.values(), dtype=float)

    ind_sort = np.argsort(form)
    ind_unsort = np.argsort(ind_sort)

    form_sort = form[ind_sort]
    grouped_ids = np.split(bhid[ind_sort], np.where(np.diff(form[ind_sort]))[0]+1)

    ngroup = 0
    now = 0
    group_sft = []
    groups = np.zeros(len(bhid),dtype=int)
    for i,g in enumerate(grouped_ids):
        if len(g)<csize:
            groups[now:now+len(g)] = ngroup + istart
            group_sft.append(form_sort[now+len(g)-1])
            now += len(g)
            ngroup +=1
        else:
            nleft = len(g)
            while nleft >= csize:
                groups[now:now+csize] = ngroup + istart
                group_sft.append(form_sort[now+csize-1])
                now += csize
                ngroup += 1
                nleft -= csize

            groups[now:now+nleft] = ngroup + istart
            group_sft.append(form_sort[now+nleft-1])
            now += nleft
            ngroup += 1
    print('Total Chunks:',ngroup,flush=True)
    return bhid[ind_sort], form_sort, groups


def append_data(bhid,sft,groups):
    path = '/hildafs/datasets/Asterix/BH_details_dict/Read-Blackhole-Detail'
    # os.makedirs(path)
    dest     = BigFile(path)
    BHID     = dest['BHID'][:]
    Index    = dest['Index'][:]
    SeedTime = dest['SeedTime'][:]
    
    
    print('Making sure we do not append twice...', flush=True)
    mask = Index < cbeg
    BHID = BHID[mask]
    Index = Index[mask]
    SeedTime = SeedTime[mask]
    
    print('N BHs before appending:',len(BHID),flush=True)
    print('Maximum seeding time before appending', max(SeedTime), flush=True)
    print('Maxumum chunk number before appending', max(Index), flush=True)
    
    
    BHID = np.concatenate([BHID,bhid])
    Index = np.concatenate([Index,groups])
    SeedTime = np.concatenate([SeedTime,sft])
    
    print('N BHs after appending:',len(BHID),flush=True)
    print('Maximum seeding time after appending', max(SeedTime), flush=True)
    print('Maxumum chunk number after appending', max(Index), flush=True)
    
    print('Writing data...',flush=True)

    blockname='Index'
    dest.create_from_array(blockname,Index)

    blockname='BHID'
    dest.create_from_array(blockname,BHID)

    blockname='SeedTime'
    dest.create_from_array(blockname,SeedTime)
    
    print('Finished writing',flush=True)

    
if __name__ == "__main__":   
    #------------------- Arguments -------------------------
    parser = argparse.ArgumentParser(description='assign file numbers for each BH based on seeding time')
    parser.add_argument('--sprev',required=True,type=int,help='last snap of the previous run')
    parser.add_argument('--sbeg',required=True,type=int,help='beginning snapshot for resuming the assignment')
    parser.add_argument('--csize',default=5000,type=int,help='number of BHs per chunk')
    parser.add_argument('--cbeg',required=True,type=int,help='first chunk number to resume the assignment')

    args = parser.parse_args()
    sbeg = int(args.sbeg)
    sprev = int(args.sprev)
    csize = int(args.csize)
    cbeg = int(args.cbeg)
    
    
    bh_sft = load_bh_sft(sprev,sbeg)
    bhid,sft,groups = assign_chunks(cbeg,csize,bh_sft)
    append_data(bhid,sft,groups)
    
    print('Done!')
