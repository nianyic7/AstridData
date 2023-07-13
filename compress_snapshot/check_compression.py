import numpy as np
from bigfile import BigFile
from bigfile import FileMPI
import sys,os
from mpi4py import MPI
import argparse
import time
import glob
import matplotlib.pyplot as plt

# we do not touch the integer columns

TypeMap = { np.dtype('float32') : np.dtype('float16'), 
            np.dtype('float64') : np.dtype('float32'), 
            np.dtype('uint64')  : np.dtype('uint64'), 
            np.dtype('uint32')  : np.dtype('uint32'), 
            np.dtype('uint8')   : np.dtype('uint8'), 
          np.dtype('int64')  : np.dtype('int64'), 
          np.dtype('int32')  : np.dtype('int32'), 
          np.dtype('int8')  : np.dtype('int8'), }

IntTypes = set([np.dtype('uint64'), np.dtype('uint32'), np.dtype('uint8'), np.dtype('int64'), np.dtype('int32'), np.dtype('int8')])

KeepCols = set(["ElectronAbundance", "Metallicity", "Metals"])

def compress(col, idata):
    if idata.dtype in IntTypes:
        return idata
    elif col == "Position":
        # direct compression
        odata = idata.astype(TypeMap[idata.dtype])
        return odata
    elif col in KeepCols: # these columns may have overflow issue during compression
        return idata
    else:
        odata = np.power(idata, 0.2).astype(TypeMap[idata.dtype])
    return odata


def recover(col, odata):
    if odata.dtype in IntTypes:
        return odata
    elif col == "Position":
        # direct compression
        rdata = odata.astype(np.dtype('float64'))
        return rdata
    elif col in KeepCols: # these columns may have overflow issue during compression
        return odata
    else:
        rdata = odata.astype(np.dtype('float32'))
        rdata = np.power(rdata, 5)
    return rdata
    


def bin_chunk(bf_r, bf_w, rank, Ntot, size, blockname, err_bins):
    beg = Ntot * rank // size
    end = Ntot * (rank  + 1) // size
    
    inner_size = nloops  
    Nthis = end - beg
    col = blockname.split('/')[-1]
    
    
    tot_counts, _ = np.histogram([], bins=err_bins)
        
    for i in range(inner_size):
        ibeg = beg + Nthis * i // inner_size
        iend = beg + Nthis * (i + 1) // inner_size
        
        idata = bf_r[blockname][ibeg : iend]
        odata = bf_w[blockname][ibeg : iend]
        rdata = recover(col, odata)

        # make sure that zeros are aligned and ignore them in the histogram
        mask = idata==0
        assert (rdata[mask] == 0).all(), "Block %s original data is 0, compressed data is not!"%blockname
        
        # check for NaNs
        if np.isnan(rdata).any():
            raise ValueError("Block %s contains NaN values after compression"%blockname)
        
        # bin relative error
        mask = idata != 0
        if not np.any(mask):
            if rank == 0:
                print('Block %s: All entries are zero'%blockname, flush=True)
            continue
            
        rel_err = (idata[mask] - rdata[mask]) / np.abs(idata[mask])
        counts, _ = np.histogram(rel_err, bins=err_bins)
        tot_counts += counts

        if i == inner_size:
            assert (iend == end)  
    return tot_counts

def save_plots(err_bins, tot_counts, blockname, vmin, vmax):
    column = blockname.split('/')[-1]
    ptype  = blockname.split('/')[-2]

    fig, ax = plt.subplots(1,1,figsize=(7,5))
    ax.bar(err_bins[:-1], tot_counts/Ntot, width=np.diff(err_bins), edgecolor="white", align="edge")
    ax.set(yscale='log', xlabel='relative error', ylabel='count', xlim=(vmin, vmax))
    ax.set_xscale('symlog', linthresh=1e-6)
    
    title = '%s-%s'%(ptype, column)
    ax.set_title(title)
    plt.savefig('%s.png'%(title), bbox_inches='tight')
    print('done saving stats for block:', blockname, flush=True)
    

if __name__ == "__main__":
    
    #----------- MPI Init ----------------
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    #-------------- Cmd line Args ------------------------------
    parser = argparse.ArgumentParser(description='PART file compression')
    parser.add_argument('--ifile',required=True,type=str,help='path of the input file (end without /)')
    parser.add_argument('--ofile',required=True,type=str,help='path of the output file directory (end without /)')
    parser.add_argument('--nloops',required=True,type=int,help='number of loops per task')
    parser.add_argument('--blocknames', nargs='+', help='blocks to process', default=None)
    
    args = parser.parse_args()
    
    #--------------------------

    
    bf_r  = BigFile(args.ifile)
    bf_w  = BigFile(args.ofile)

    nloops = int(args.nloops)
    
    # block names to process
    blocknames = []
    if args.blocknames is None:
        all_blocks = sorted(glob.glob(args.ifile + '/*/*'))
        for ff in all_blocks:
            column = ff.split('/')[-1]
            ptype  = ff.split('/')[-2]
            blockname = ptype + '/' + column

            if (ptype in ['0', '1', '4', '5']):
                blocknames.append(blockname)
    else:
        for bn in args.blocknames:
            if os.path.isdir(args.ifile + '/' + bn):
                blocknames.append(bn)
            elif rank == 0:
                print('Block %s does not exit!'%(bn), flush=True)

    comm.barrier()
    
    if rank == 0:
        print('All blocks to check:', flush=True)
        for bb in blocknames:
            print(bb, flush=True)
    
    comm.barrier()
    if rank == 0:
        t1 = time.time()
        
    # create figure for plotting
    nblocks = len(blocknames)
    # nrows = nblocks//5+1
    # fig, ax = plt.subplots(nrows, 5, figsize=(nrows*5, 5*5), sharex=True)
    # ax = ax.flatten()
    
    pos_bins = np.logspace(-6,0.2,20)
    neg_bins = -pos_bins[::-1]
    
    err_bins = np.concatenate([neg_bins, pos_bins])
    print(err_bins)
    vmin, vmax = min(err_bins), max(err_bins)
    
    for blockname in blocknames:
        Ntot = bf_r[blockname].size
        tot_counts = bin_chunk(bf_r, bf_w, rank, Ntot, size, blockname, err_bins)
        comm.barrier()
        
        tot_counts = comm.allreduce(tot_counts,op=MPI.SUM)

        if rank == 0:
            save_plots(err_bins, tot_counts, blockname, vmin, vmax)
            t2 = time.time()
            print("Finished checking {: <15}".format(blockname), flush=True)
            print("time used = %.1f min:"%((t2 - t1)/60.), flush=True)
            t1 = t2	
        comm.barrier()
    
    if rank == 0:
        print("Done!", flush=True)
