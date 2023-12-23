import numpy as np
import h5py
import glob
import argparse
from tempfile import mkstemp
from shutil import move, copymode
from os import fdopen, remove
import os
import pickle
import numpy.lib.recfunctions as rf 

root = '/jet/home/nianyic/scratch1/asterix_binary/'

binaries = np.load(root+'binary_parameters_vcut=5.0_kpc_rextr=1.2_soft_wseed.npy')
print('Loaded binaries, N=%d'%len(binaries))



root = '/jet/home/nianyic/scratch1/asterix_binary/binary_galaxy-z3/'

dpre = {}
for subdir in sorted(glob.glob(root+'*tomerge*.pkl')):
    print('Processing:',subdir)
    with open(subdir, 'rb') as f:
        subdict_pre = pickle.load(f)
    f.close()
    dpre.update(subdict_pre)
print('Loaded all binaries to merge:',len(dpre))

dpost = {}
for subdir in sorted(glob.glob(root+'*merged*.pkl')):
    print('Processing:',subdir)
    with open(subdir, 'rb') as f:
        subdict_post = pickle.load(f)
    f.close()
    dpost.update(subdict_post)
print('Loaded all merged binaries:',len(dpost))

      

# update binary variables
nmiss = 0
f2 = ['id1','id2','zmerge','mass1','mass2','seed1','seed2','sig1','sig2','sigma','rho','gamma','mstot','reff']
dtype2 = ['q','q','d','d','d','d','d','d','d','d','d','d','d','d']
dt1 = {'names':f2, 'formats':dtype2}
dt = np.dtype(dt1)
new_data = []
for b in binaries:
    id1 = b['id1']
    id2 = b['id2']
    key = (id1,id2)
    try:
        bpre = dpre[key]
    except KeyError:
        print('binary not found in previous snap:',key)
        nmiss += 1
        continue
    try:
        bpost = dpost[key]
    except KeyError:
        print('binary not found in next snap:',key)
        nmiss += 1
        continue

    data = [b['id1'],b['id2'],b['zmerge'],b['mass1'],b['mass2'],b['seed1'],b['seed2'],\
            bpre['sigma1'],bpre['sigma2'],bpost['sigma'],bpost['rho'],bpost['gamma'],\
            bpost['mstot'],bpost['reff']]

    new_data.append(data)
    
new_data = rf.unstructured_to_structured(np.array(new_data), dt) 
print('N missed binaries:',nmiss)      
print('data shape:',new_data.shape)  
print('data,dtype:',new_data.dtype)
print(new_data[:3])
print(new_data['id1'][:3])
print('Saving all binary info...')

root = '/jet/home/nianyic/scratch1/asterix_binary/'      
r = 'binary_subfind-z3'
np.save(root+r,new_data)
      
     
r = 'binary_info_pre-z3'
with open(root+r+'.pkl', 'wb') as f:
    pickle.dump(dpre,f)
    f.close()   
         
r = 'binary_info_post-z3'
with open(root+r+'.pkl', 'wb') as f:
    pickle.dump(dpost,f)
    f.close()             
      
         
