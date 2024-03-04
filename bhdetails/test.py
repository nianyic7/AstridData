import os, glob,struct
import numpy as np
fields = ('size','BHID','BHMass','Mdot','Density','timebin','encounter','MinPos',\
         'MinPot','Entropy','GasVel','acMom','acMass','acBHMass',\
         'Fdbk','SPHID','SwallowID','CountProgs','Swallowed',\
         'BHpos','srDensity','srParticles','srVel','srDisp',\
         'DFAccel','DragAccel','GravAccel','BHvel','Mtrack','Mdyn',\
         'KinFdbk','NumDM','V1sumDM','V2SumDM','MgasEnc','KEflag','z','size2')
dtype = ('i','q','d','d','d','i','i','3d',\
        'd','d','3d','3d','d','d',\
        'd','q','q','i','i',\
        '3d','d','d','3d','d',\
        '3d','3d','3d','3d','d','d',\
        'd','d','3d','d','d','i','d','i')
dt = {'names':fields, 'formats':dtype}

data = np.zeros(1, dtype=dt)
nsingle = data.nbytes


def read_binary_file(file):
    fields = ('size','BHID','BHMass','Mdot','Density','timebin','encounter','MinPos',\
             'MinPot','Entropy','GasVel','acMom','acMass','acBHMass',\
             'Fdbk','SPHID','SwallowID','CountProgs','Swallowed',\
             'BHpos','srDensity','srParticles','srVel','srDisp',\
             'DFAccel','DragAccel','GravAccel','BHvel','Mtrack','Mdyn',\
             'KinFdbk','NumDM','V1sumDM','V2SumDM','MgasEnc','KEflag','z','size2')
    dtype = ('i','q','d','d','d','i','i','3d',\
            'd','d','3d','3d','d','d',\
            'd','q','q','i','i',\
            '3d','d','d','3d','d',\
            '3d','3d','3d','3d','d','d',\
            'd','d','3d','d','d','i','d','i')
    dt = {'names':fields, 'formats':dtype}
    #-----------------------------------------------------------
    struct_fmt = ''.join(list(dtype))
    struct_len = struct.calcsize(struct_fmt)
    data = np.fromfile(file, dtype=dt, count=-1)
    results = data

    results['z'] = 1./results['z'] - 1
    print("read %d BHs"%(len(results)))
    print(f"Total memory occupied: {results.nbytes/1e9} GBs")
    return results

file_name = "/hildafs/datasets/Asterix/BH_details_bigfile2/test"
flist = sorted(glob.glob(file_name + "/*"))

for ff in flist[:20]:
    file_stats = os.stat(ff)
    print(ff)
    print(file_stats)
    print(f"File Size in Bytes is {file_stats.st_size}")
    print(f"number of element {file_stats.st_size / nsingle}")
    read_binary_file(ff)
