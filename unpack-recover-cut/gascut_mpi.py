import numpy as np
from bigfile import BigFile
import os, sys
from mpi4py import MPI


center = np.array([80894.68825668, 820.77664037, 71608.47575295]) # position of BH0
BoxSize = 250000.0 # kpc/h
crop = 10000.0 # kpc/h

def main():
    
    part = BigFile('/home1/08942/nianyic/scratch3/Astrid/PART_temp/PART_021')
    
    comm = MPI.COMM_WORLD
    
    prefix = '/home1/08942/nianyic/scratch3/Astrid/PART_temp/cut/'
    ofilename = prefix + 'PART_021_4'
    
    if comm.rank == 0:
        ofile = BigFile(ofilename, create=True)
        print('writing to ', ofilename)
        
    comm.barrier()
    
    ofile = BigFile(ofilename, create=True)
    
    #######################################################
    Ntot = 5500**3
    
    start = Ntot * comm.rank // comm.size
    end = Ntot * (comm.rank  + 1) // comm.size     
    
    Nchunks = 800 # split the particles into Nchunks of groups for reading position    
    chunksize = int(Ntot/Nchunks) # each group has chunksize of particles
    
    gaspos = []
    gasm = []
    u = []
    ne = []
    sml = []
    nhf = []
    
    for i in range(start,end,chunksize):
        sl = slice(i,i+chunksize)
        
        pos = np.float32(part['0/Position'][sl])
        pos -= center

        pos[pos < - BoxSize*0.5] += BoxSize
        pos[pos > BoxSize*0.5] -= BoxSize
        
        mask = np.abs(pos[:,0])  < crop
        mask &= np.abs(pos[:,1]) < crop
        mask &= np.abs(pos[:,2]) < crop
        
        ld = mask.sum()
        
        if ld>0:
            gaspos.append(part['0/Position'][sl][mask])
            u.append(part['0/InternalEnergy'][sl][mask])
            ne.append(part['0/ElectronAbundance'][sl][mask])
            sml.append(part['0/SmoothingLength'][sl][mask])
            gasm.append(part['0/Mass'][sl][mask])
            nhf.append(part['0/NeutralHydrogenFraction'][sl][mask])
    
    if len(gasm)>0:
        gaspos = np.array(np.concatenate(gaspos))
        gasm = np.array(np.concatenate(gasm))
        u = np.array(np.concatenate(u))
        ne = np.array(np.concatenate(ne))
        sml = np.array(np.concatenate(sml))   
        nhf = np.array(np.concatenate(nhf)) 
    
        gaspos -= center
        gaspos[gaspos < - BoxSize*0.5] += BoxSize
        gaspos[gaspos > BoxSize*0.5] -= BoxSize
        gaspos += np.array([crop,crop,crop])
    
    comm.barrier()
    ###############################################
    olength = len(gasm)
    ooffset = sum(comm.allgather(olength)[:comm.rank])
    oN = comm.allreduce(olength)
    
    blocks = ['0/Position', '0/Mass', '0/InternalEnergy',
              '0/ElectronAbundance','0/SmoothingLength',
              '0/NeutralHydrogenFraction']
    data = [gaspos,gasm,u,ne,sml,nhf]
    
    for i,block in enumerate(blocks):
        iblock = part[block]
        if comm.rank == 0:
            print(block, 'size', oN)
            ofile.create(block, size=oN, dtype=iblock.dtype, Nfile=8)
        comm.barrier()
        oblock = ofile[block]
        oblock.write(ooffset, np.array(data[i]))
        

if __name__ == "__main__":
    main()
