from bigfile import BigFile
import numpy as np
import sys
import h5py
import sys,os
import pickle
import argparse

class NameMaps:
    """Maps names between HDF5 and Bigfile."""
    def __init__(self):
        #Map of (only changed) block names between HDF5 and bigfile snapshots.
        self.hdf_to_bigfile_map = { "Coordinates" : "Position",
                                    "Velocities": "Velocity", 
                                    "Masses": "Mass",
                                    "NeutralHydrogenAbundance": "NeutralHydrogenFraction",
                                    "GFM_Metallicity": "Metallicity",
                                    "ParticleIDs" : "ID",
                                  }
        self.bigfile_to_hdf_map = {v : k for (k,v) in self.hdf_to_bigfile_map.items()}
        #Leave the metallicity array unchanged.
        del self.bigfile_to_hdf_map["Metallicity"]

    def get_bigfile_name(self, hdfname):
        """Get the bigfile name from an HDF5 name."""
        try:
            return self.hdf_to_bigfile_map[hdfname]
        except KeyError:
            #Unrecognised names are assumed unchanged
            return hdfname

    def get_hdf5_name(self, bigfilename):
        """Get the HDF5 name from a bigfile name."""
        try:
            return self.bigfile_to_hdf_map[bigfilename]
        except KeyError:
            #Unrecognised names are assumed unchanged
            return bigfilename

#Global name registry
names = NameMaps()

bf_blocks = ['0/ID','0/Position','0/Velocity','0/Mass','0/InternalEnergy',
             '1/ID','1/Position','1/Velocity','1/Mass',
             '4/ID','4/Position','4/Velocity','4/Mass',
             '5/ID','5/Position','5/Velocity','5/Mass',]

def compute_nfiles(npart):
    """Work out how many files we need to split the snapshot into.
       We want less than 2^31 bytes per data array, and we want a power of two,
       so that we probably divide the particles evenly."""
    nfiles = 1
    #Largest possible data array: a 3-vector in double precision.
    maxarray = np.max(npart) * 3 * 8
    while maxarray // nfiles >= 2**31:
        nfiles *=2
    return nfiles

def write_hdf_header(bf, hdf5, npart_tot, nfiles, npart_file):
    """Generate an HDF5 header from a bigfile PIG header."""
    head = hdf5.create_group("Header")
    hattr = head.attrs
    battr = bf["Header"].attrs
    #As a relic from Gadget-1, the total particle numbers
    #are written as two separate 32 bit integers.
    hattr["NumPart_Total"] = np.uint32(npart_tot % 2**32)
    hattr["NumPart_Total_HighWord"] = np.uint32(npart_tot // 2**32)
    hattr["NumPart_ThisFile"] = np.int32(npart_file)
    hattr["NumFilesPerSnapshot"] = np.int32(nfiles)
    hattr["Redshift"] = 1./battr["Time"] - 1
    hattr["MassTable"] = np.int32([0,0,0,0,0,0]) # we input mass instead
    #Pass other keys through unchanged. We whitelist expected keys to avoid confusing Gadget.
    hdfats = ["Time", "BoxSize", "Omega0", "OmegaLambda", "HubbleParam", "OmegaBaryon"]
    for attr in hdfats:
        hattr[attr] = battr[attr]

def write_hdf_file(bf, hdf5name,groups,nfiles):
    """Write the data arrays to an HDF5 file."""
    ostart = obt[groups]
    oend = obt[groups+1]
    npart = np.sum(lbt[groups],axis=0) # NumPartInGroupTotal
    # Variables for the velocity factors
    atime = bf["Header"].attrs["Time"]
    perfile = (npart//nfiles)
    print('Particles per file:',perfile)
    
    for fnum in range(nfiles):
        if nfiles==1:
                hfile = hdf5name +".hdf5"
        else:
            hfile = hdf5name + "."+str(n)+".hdf5"
        with h5py.File(hfile,'w') as hdf5:      
            #Write the header
            write_hdf_header(bf, hdf5, npart, nfiles, npart//nfiles)       
            for ptype in [0,1,4,5]:
                hdf5.create_group("PartType" + str(ptype))
                
    for block in bf_blocks:
        ptype, bname = os.path.split(block)
        ptype = int(ptype)
        nn = perfile[ptype]
        hname = names.get_hdf5_name(bname)
        # Load all needed data in this field
        bfdata = np.concatenate([bf[block][ostart[i][ptype]:oend[i][ptype]] for i in range(len(ostart))],axis=0)
        print(block)
        print('data shape:',bfdata.shape,flush=True)
        if (bname == "Velocity"):
            bfdata /= np.sqrt(atime) #pecvel is a dx/dt: (physical peculiar velocity)   
        if (bname == "InternalEnergy"):
            bfdata = np.zeros(npart[0])
        # Now split into nfiles
        for n in range(nfiles):
            if nfiles==1:
                hfile = hdf5name +".hdf5"
            else:
                hfile = hdf5name + "."+str(n)+".hdf5"
            print('saving from %d to %d'%(n*nn,(n+1)*nn))
            if n == nfiles - 1:
                bb = bfdata[int(n*nn):]
            else:
                bb = bfdata[int(n*nn):int((n+1)*nn)]
            with h5py.File(hfile,'a') as hdf5:      
                #Write the data
                hdf5["PartType"+str(ptype)][hname] = bb
        print("Wrote file %d" % n, 'to '+hfile)


def get_obt(pig):
    lbt = pig.open('FOFGroups/LengthByType')[:]
    obt = np.cumsum(lbt,axis=0)
    a1 = np.array([[0,0,0,0,0,0]],dtype=np.uint64)
    obt = np.append(a1,obt,axis=0)
    return obt

def split_group(pig,groupids,nmax=512):
    """
    Find BHIDs of binaries that will merge in the next snapshot
    i.e. this is the last snapshot before the merger
    snapz[i] > zmerge >= snapz[i+1]
    we save both IDs since the merger has not happened

    Args:
        mergers (ndarray):        File(s) that stores the raw simulation mergers
        snaps (ndarray):          (nsnaps,2) array of all available snaps and their redshifts
        i (int):                  index of the snapshot in the array of snaps
    Returns:
        ids (set):                set of bhids to be merged soon after this snapshot
    """
    
    
    Nps = np.sum(lbt[groupids-1],axis=0) # NumPartInGroupTotal
    TNs = np.sum(Nps)

    Njob = np.int(TNs/(nmax**3))+1
    Nchunk = int(TNs/Njob)+1


    lbt_all = np.sum(lbt[groupids-1],axis=-1)
    off_all = np.cumsum(lbt_all)

    # off_all = np.append(np.array([0]),off_all)
    idxs = []
    for j in range(0,Njob):
        ofs = j*Nchunk
        mask = off_all>ofs
        mask &= off_all<= ofs+Nchunk
        gind = groupids[mask]-1
        idxs.append(gind)
    return idxs

def get_pig(snap):
    if snap > 294:
        pig_dir = astrid_path + 'PIG2/PIG_%03d'%snap
        assert os.path.isdir(pig_dir), "%s does not exist!"%pig_dir
        pig = BigFile(pig_dir)
    else:
        pig_dir = astrid_path + 'PIG_files/PIG_%03d'%snap
        assert os.path.isdir(pig_dir), "%s does not exist!"%pig_dir
        pig = BigFile(pig_dir)
    return pig




if __name__ == "__main__":
    
    # GLOBAL INFO
    astrid_path = '/hildafs/datasets/Asterix/'
    group_path  = '/hildafs/home/nianyic/scratch1/Astrid_data/binary_cat/subfind/grpList/'
    pig_path    = astrid_path + 'PIG*/'

    #------------------- Arguments --------------------------------------
    parser = argparse.ArgumentParser(description='subgrpID-subfind')
    parser.add_argument('--snap',required=True,type=int,help='path of the subfind directory')
    parser.add_argument('--savedir',required=True,type=str,help='directory to save hdf files')

    args = parser.parse_args()
    snap = int(args.snap)
    savedir = args.savedir

    #----------------------------------------------------------------------------------#

    pig = get_pig(snap)
    # get obt
    lbt = pig['FOFGroups/LengthByType'][:]
    obt = get_obt(pig)

    groupids = np.load(group_path + 'merger-group-%03d.npy'%snap)
    groups = split_group(pig,groupids,nmax=512)

    # Output some basic group info
    print ('Groups in each chunk:', [len(k) for k in groups])
    print('Total groups:',sum(np.array([len(k) for k in groups])))
    print('Total chunks:', len(groups))
    

    for j,g in enumerate(groups):
        print('Processing chunk:',j)
        ostart = obt[g]
        oend = obt[g+1]
        npart = np.sum(lbt[g],axis=0) # NumPartInGroupTotal
        print('Particles in this chunk:',npart)
        nfiles = compute_nfiles(npart)

        print("Writing %d hdf snapshot files" % (nfiles),flush=True)

        outputdir = savedir+'snap%03d'%snap+'/chunk{}.{}/output'.format(j,g[-1])
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)
        ofilename = outputdir+'/snap_%03d'%snap
        write_hdf_file(pig, ofilename,g,nfiles)

    