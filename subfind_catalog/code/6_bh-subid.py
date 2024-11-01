import numpy as np
from bigfile import BigFile,FileMPI
import sys,os,gc
import glob
import argparse


# example to run 
# python3 $codedir/6_bh-subid_v3.py  --dest $indir --particle_type 5 


if __name__ == "__main__":
    

    #-------------- Cmd line Args ------------------------------
    parser = argparse.ArgumentParser(description='subgrpID-subfind')
    parser.add_argument('--src',required=True,type=str,help='path of the input file directory')
    parser.add_argument('--dest',required=True,type=str,help='path of the output file directory')
    parser.add_argument('--particle_type',required=True,type=int,help='the particle type to be processed')
    args = parser.parse_args()
    
    pig  = BigFile(args.src)
    pig_w = BigFile(args.dest, create=True)
    particle_type = int(args.particle_type)    


    
    #-------------- Initialize the block ------------------------------
    dtype = 'i8'
    dsize  = pig[f'{particle_type}/GroupID'].size
    nfile = pig[f'{particle_type}/GroupID'].Nfile

    blockname = f'{particle_type}/SubgroupIndex'
    hostgal_block = pig_w.create(blockname, dtype, dsize, nfile)
    print('Initizlized block:', blockname,flush=True)
    print('size:',dsize, 'type:',dtype,flush=True)    

    

    suboff = pig["SubGroups/SubhaloOffsetType"][:][:,particle_type]
    sublen = pig["SubGroups/SubhaloLenType"][:][:,particle_type]
    part_num = pig["Header"].attrs["NumPartInGroupTotal"][particle_type]
    # last_particle is the last particle included in the subgroup
    last_particle = suboff[-1] + sublen[-1]
    print(f"{particle_type} Total number of particles: {part_num} last_particle {last_particle}", flush=True)
    
    
    istart = 0
    chunk_size = int(1e9)
    FUZZ = 100000000
    print(f"chunksize {chunk_size}", flush=True)

    while(istart < last_particle):
        iend = istart + chunk_size
        iend = min(iend, last_particle)
        print(particle_type, istart, iend, "Start", flush=True)
        bh_idx_list = np.arange(istart,iend)
        subhalo_idx = np.searchsorted(suboff,bh_idx_list, side='right')-1
        
        #-------------- deal with the "fuzz" particle ------------------------------
        # fuzz particles were labeled in SubgroupID (written in post-process step1)
        # in the new Subgroup_Index block, fuzz particle are labeled as "-1"
        subhalo_chunkidx = pig[f"{particle_type}/SubgroupID"][istart:iend]
        fuzz_mask = subhalo_chunkidx == FUZZ
        subhalo_idx[fuzz_mask] = -1
        subhalo_idx.astype(np.uint32)
        hostgal_block.write(istart, subhalo_idx)
        istart = iend
        
        del bh_idx_list, subhalo_idx, subhalo_chunkidx, fuzz_mask
        gc.collect()
        
        
    # fill the last fuzz particle 
    subhalo_idx = np.ones(int(part_num) - int(last_particle), dtype=np.uint32) * -1
    hostgal_block.write(iend, subhalo_idx)

    print(f"Particle {particle_type} Done", flush=True)
    
    
    
    
    
