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
parser.add_argument('--srcdir',default="/hildafs/datasets/Asterix/BH_details_bigfile",type=str,help='directory name of the bhdetails in bigfile')
parser.add_argument('--outdir',default="/hildafs/datasets/Asterix/BH_details_dict/first_pass",type=str,help='output dir of the bhdetails as dicts')

args = parser.parse_args()
idx = int(args.idx)
srcdir = args.srcdir
outdir = args.outdir

#--------------------------------------------------------------------------------
    
name_all = ('size','BHID','BHMass','Mdot','Density','timebin','Encounter','MinPos',\
        'MinPot','Entropy','GasVel','acMom','acMass','acBHMass',\
        'Fdbk','SPHID','SwallowID','CountProgs','Swallowed',\
        'BHpos','srDensity','srParticles','srVel','srDisp',\
        'DFAccel','DragAccel','GravAccel','BHvel','Mtrack','Mdyn',\
        'KineticFdbkEnergy','NumDM','V1sumDM','V2sumDM','MgasEnc','KEflag','z','size2')
dtype_all = ('i','q','d','d','d','i','i','3d',\
    'd','d','3d','3d','d','d',\
    'd','q','q','i','i',\
    '3d','d','d','3d','d',\
    '3d','3d','3d','3d','d','d',\
    'd','d','3d','d','d','i','d','i')
name2type = {name: dtype for name, dtype in zip(name_all, dtype_all)}

features = ["acBHMass", "BHMass", "BHpos", "BHID", "BHvel", "Mdot", "SwallowID", "z"]
dtype = [name2type[f] for f in features]

#----------------------------------------------------------------------------------
start = time.time()
t1 = time.time()

flist = sorted(glob.glob(srcdir+'/BH-Details-R%03d*'%idx))
for details in flist:
    print('Processing: ', details,flush=True)
    fname = os.path.basename(details).split('.')[0]
    bf = BigFile(details)
    bhids = bf.open('BHID')[:]
    entry = len(bhids)
    print('Total data length: %e'%entry)

    all_info = np.zeros(len(bhids),dtype=[(f,dtype[i]) for i,f in enumerate(features)])

    # read in all data
    for ff in features:
        print('Reading:',ff,flush=True)
        print('time elapsed:',time.time()-t1,flush=True)
        t1 = time.time()
        try:
            size = bf[ff].size
            half = size//2
            all_info[ff][:half] = bf[ff][:half]
            all_info[ff][half:] = bf[ff][half:]
    #        all_info[ff] = bf.open(ff)[:]
        except:
            print('field:'+ff+' does not exist!')
            continue

    print('Sorting...',flush=True)
    all_info.sort(order=["BHID", "z"])
    print('time elapsed:',time.time()-t1,flush=True)
    t1 = time.time()

    print('Splitting...',flush=True)
    all_info = np.split(all_info, np.where(np.diff(all_info['BHID']))[0]+1)
    print('total BHs:',len(all_info),flush=True)
    print('time elapsed:',time.time()-t1,flush=True)

    t1 = time.time()

    print('Creating dictionary...',flush=True)
    all_info = {b['BHID'][0]:b for b in all_info}
    print('time elapsed:',time.time()-t1,flush=True)
    t1 = time.time()

    # save
    save = outdir+f'/{fname}.pkl'
    print('Saving to:', save, flush=True)

    with open(save, 'wb') as f:
        pickle.dump(all_info,f)
        f.close()
        
    print('Done! Total time elapsed:',time.time() - start)




    
    
