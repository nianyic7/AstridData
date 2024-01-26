import numpy as np
from bigfile import BigFile
from bigfile import FileMPI
import sys,os
from mpi4py import MPI


def write_arr(rank,start,length):
    data = rank*np.ones((length,2))
    dest[blockname].write(start,data)
    # block.write(start,data)
    print('Finished writing:',start,flush=True)

if __name__ == "__main__":
    
    #----------- MPI Init ----------------
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    
    loop = 100
    length = 100000
    
    fname = '/hildafs/home/nianyic/scratch1/test_bigfile'
    dest = FileMPI(comm, fname, create=True)
    
    print('created file:',fname,flush=True)

    dtype = ('f4', (2,))
    nfile = 16
    dsize = size * length * loop
    blockname = '0/test2'
    
    
    # dest.create_from_array('0/test1', np.zeros((dsize,2)))


    comm.barrier()
    block = dest.create(blockname, dtype, dsize, nfile)

    comm.barrier()
    
    offset = rank * length * loop

    for i in range(80,loop):
        print('rank %d, iteration %d, start'%(rank,i),flush=True)
        start = offset + i*length

        write_arr(rank,start,length)
        print('rank %d, iteration %d, end'%(rank,i),flush=True)
        
        
        comm.barrier()
        
        if rank == 0:
            print('-------------- iteration done ------------, %d'%i, flush=True)
            
        comm.barrier()