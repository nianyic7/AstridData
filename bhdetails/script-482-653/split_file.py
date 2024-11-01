import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from bigfile import BigFile
import glob,os,struct
import pickle
import argparse
import time


#------------------- Arguments -------------------------
parser = argparse.ArgumentParser(description='save bhdetails as dictionaries')

parser.add_argument('--idx',required=True,type=int,help='index of the detail file')

args = parser.parse_args()
idx = int(args.idx)


#----------------------- Load dictionary -----------------------------------------
start = time.time()
t1 = time.time()

first_dir = '/hildafs/datasets/Asterix/BH_details_dict/first_pass/'
detail = first_dir+'BlackholeDetails_%03d.pkl'%idx
print('Loading from:',detail,flush=True)

with open(detail, 'rb') as f:
    data = pickle.load(f)
f.close()

print('Finished loading data',time.time() - t1,flush=True)
t1 = time.time()

#------------ figure out the starting point ------------------------
root = '/hildafs/datasets/Asterix/'
snapz_all = np.loadtxt(root+'PIG2/Snapshots.txt',unpack=True)
snapz_all[1] = 1./snapz_all[1]-1
snapz = []

for i in snapz_all[0]:
    if os.path.isfile(root+'BH_details_dict/first_pass/BlackholeDetails_%03d.pkl'%i):
        snapz.append([i,snapz_all[1,int(i)]])
snapz = np.array(snapz).T

snext = (snapz[0]==idx).nonzero()[0][0]+1
if snext==len(snapz[0]):
    zmin = 1.0    
else:  
    zmin = snapz[1][snext]
print('Minimum z in this file:',zmin,flush=True)

#--------------------- load BH groups -----------------------------------
path = '/hildafs/datasets/Asterix/BH_details_dict/Read-Blackhole-Detail'
pig = BigFile(path)
BHIDs = pig.open('BHID')[:]
Index = pig.open('Index')[:]
ST = pig.open('SeedTime')[:]

splits = np.where(np.diff(Index))[0]+1
group_id = np.split(BHIDs, splits)
group_seed = np.concatenate([ST[0:1],ST[splits]])
print('number of chunks to split into:', len(group_id), flush=True)


#--------------------- Loop over groups ---------------------------------
outdir = '/hildafs/datasets/Asterix/BH_details_dict/second_pass/'

for chunk,bhid in enumerate(group_id):
    zseed = 1./group_seed[chunk] - 1.
    
    # skip BHs seeded after this redshift
    if zseed < zmin:
        print('skipping chunk:', chunk, 'because zseed < zmin', flush=True)
        continue
    print('Chunk:',chunk,'Number of BHs:',len(bhid),flush=True)
    print('Seeding time for this chunk:',zseed,flush=True)
    bhs = {}
    for b in bhid:
        try:
            bhs[b] = data[b]
        except KeyError:
            continue
            
    # save
    save = outdir+'BlackholeDetails-%03d-%04d.pkl'%(idx,chunk)
    print('Saving to:', save, flush=True)

    with open(save, 'wb') as f:
        pickle.dump(bhs,f,pickle.HIGHEST_PROTOCOL)
    f.close()
    
    print('Finished chunk ', chunk, time.time()-t1,flush=True)
    t1 = time.time()

print('Done! Total time:',time.time() - start,flush=True)


    
    
