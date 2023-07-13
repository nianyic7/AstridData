import bigfile as bf
import numpy as np
from astropy import cosmology, units as U, cosmology as C

from nbodykit.lab import BigFileCatalog
import dask.array as da
from mpi4py import MPI
from matplotlib import pyplot as plt

import nbodykit
nbodykit.use_mpi()
nbodykit.set_options(dask_chunk_size=1024*1024*8)
nbodykit.setup_logging(log_level='info')

#nbodykit has annoying deprecation warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import logging
logger = logging
logging.basicConfig(level=logging.INFO)

cosmo = cosmology.FlatLambdaCDM(Om0=0.2814,Ob0=0.0464,H0=69.7)
mean_bary_dens_cgs = (cosmo.critical_density(0) * cosmo.h**(-2) * cosmo.Ob0).to('g cm-3').value

#/*convert U (erg/g) to T (K) : U = N k T / (γ - 1)
#T = U (γ-1) μ m_P / k_B
#where k_B is the Boltzmann constant
#γ is 5/3, the perfect gas constant
#m_P is the proton mass
#μ is 1 / (mean no. molecules per unit atomic weight) calculated in loop.
#Internal energy units are 10^-10 erg/g*/
def u_to_t(uin,xhi):
    helium = 0.24
    #assuming hei ion with HI
    nep = (1-3/4*helium)*(1 - xhi)
    hy_mass = 1 - helium
    muienergy = 4 / (hy_mass * (3 + 4*nep) + 1)*uin
    temp = 5/3 * 1.6726e-24 / 1.38066e-16 * muienergy * 1e10
    return temp

def get_snap_list(outdir):
    time_list = np.loadtxt(f'{outdir}/Snapshots.txt', dtype=float, ndmin=2)
    snapshot_list = time_list[:,0]
    time_list = time_list[:,1]
    redshift_list = 1/time_list - 1
    return snapshot_list,time_list,redshift_list

def parse_property(dset, prop, **kwargs):
    #special cases where we need a slice, unit changes or transformation
    if prop == 'Temperature':
        arr = u_to_t(dset['InternalEnergy'],dset['NeutralHydrogenFraction'])
    elif prop == 'Overdensity':
        arr = dset['Density'] * dset.attrs['UnitMass_in_g'] / dset.attrs['UnitLength_in_cm']**3
        arr = (arr / mean_bary_dens_cgs) #should be unitless
    elif prop == 'Pos_flat':
        arr = dset['Position'].reshape((dset["Position"].size,))
    elif prop == 'Vel_flat':
        arr = dset['Velocity'].reshape((dset["Velocity"].size,))
    else:
        try:
            arr = dset[prop]
        except:
            logger.warning(f"cannot find {prop}")
            return None

    return arr

def print_property_sums(bfile,ptype='0/'):
    dset = BigFileCatalog(f'{bfile}',dataset=ptype)
    comm = dset.comm
    if comm.rank == 0:
        logger.info(f'Printing sums for file {bfile}')
        logger.info(f'-----------------------------')
    for col in dset.columns:
        dsum = dset[col].sum().compute()
        csum = comm.allreduce(dsum,op=MPI.SUM)

        dmin = dset[col].min().compute()
        cmin = comm.allreduce(dmin,op=MPI.MIN)
        
        dmax = dset[col].max().compute()
        cmax = comm.allreduce(dmax,op=MPI.MAX)
        
        cavg = csum / dset.csize
        if comm.rank == 0:
            logger.info(f'Property "{col}" has sum {csum:.8e} avg {cavg:.8e} min {cmin:.8e} max {cmax:.8e}')


def print_property_diffs(bfile1,bfile2,ptype='0/'):
    dset1 = BigFileCatalog(f'{bfile1}',dataset=ptype)
    comm = dset1.comm
    dset2 = BigFileCatalog(f'{bfile2}',dataset=ptype)
    if comm.rank == 0:
        logger.info(f'Printing log differences for files {bfile1}, {bfile2}')
        logger.info(f'-----------------------------')
    for col in dset1.columns:
        if col in ["Weight","Value","Selection"]:
            continue

        #NOTE: ASSUMING IDENTICAL ID ORDER
        #logdiff = da.log10(dset2[col]/dset1[col])
        diff = 2*(dset2[col] - dset1[col])/(dset2[col] + dset1[col])
        
        #only nans occur due to both datasets == zero, since these are the same I want to count them as no difference
        #+- infinity shouldn't occur, since we aren't logarithming anything
        diff = da.nan_to_num(diff)
        dmean = da.mean(diff).compute()
        drms = da.sqrt(da.mean(diff**2)).compute()
        cmean = comm.allreduce(dmean,op=MPI.SUM)
        crms = comm.allreduce(drms,op=MPI.SUM)
        
        if comm.rank == 0:
            logger.info(f'Property "{col}" has mean fractional diff {cmean:.8e} rms {crms:.8e}')

#properties should be a list of strings containing the columns we want to compare
#limits should be a ndarray (n_properties,2) or something that can be turned into one
#log should be a boolean list of each property to be binned in logspace (True) or linspace (False)
def get_particle_ndhist(bfile,properties,limits,log,ptype='0/',nbins=50):
    logger.info(f'making particle ndhist for {properties}')
    try:
        limits = np.asarray(limits)
    except:
        raise ValueError("ndhist limits should be array-like")

    if limits.shape != (len(properties),2):
        raise ValueError("ndhist limits should have shape (n_properties,2)")

    dset = BigFileCatalog(f'{bfile}',dataset=ptype)
    comm = dset.comm
    edges = np.array([np.logspace(np.log10(l[0]),np.log10(l[1]),num=nbins) if o else np.linspace(l[0],l[1],num=nbins) for l,o in zip(limits,log)])

    #NOTE: testing different ways to setup the dataset 
    #data = da.stack([dset[prop] for prop in properties])
    #data = data.rechunk((100000,4))
    
    data = [parse_property(dset,prop) for prop in properties]
    buf,edges = da.histogramdd(data,edges)

    hist = buf.compute()

    #each rank has a subset of particles, so summing the histograms will give the final result
    hist = comm.allreduce(hist,op=MPI.SUM)

    return hist,edges

#same as above but returns list of 1d histograms
def get_particle_1dhist(bfile,properties,limits,log,ptype='0/',nbins=50):
    logger.info(f'making particle 1dhists for {properties}')
    try:
        limits = np.asarray(limits)
    except:
        raise ValueError("ndhist limits should be array-like")

    if limits.shape != (len(properties),2):
        raise ValueError("ndhist limits should have shape (n_properties,2)")

    dset = BigFileCatalog(f'{bfile}',dataset=ptype)
    comm = dset.comm
    edges_in = np.array([np.logspace(np.log10(l[0]),np.log10(l[1]),num=nbins) if o else np.linspace(l[0],l[1],num=nbins) for l,o in zip(limits,log)])

    hist_arr = np.zeros((len(properties),nbins-1))
    edge_arr = np.zeros((len(properties),nbins))
    for i,prop in enumerate(properties):
        logger.info(f'starting {prop}')
        data = parse_property(dset,prop)
        buf,edges = da.histogram(data,edges_in[i])

        hist = buf.compute()

        #each rank has a subset of particles, so summing the histograms will give the final result
        hist = comm.allreduce(hist,op=MPI.SUM)
        hist_arr[i,:] = hist
        edge_arr[i,:] = edges

    return hist_arr,edge_arr

#logdiff histogram between two datasets
#NOTE: UNFINISHED FOR COMPARING MODELS/SNAPSHOTS, JUST USE FOR COMPRESSION TEST
def get_diff_ndhist(bfiles,properties,limits,ptype='0/',nbins=100):
    logger.info(f'making log difference ndhist for {properties}')
    try:
        limits = np.asarray(limits)
    except:
        raise ValueError("ndhist limits should be array-like")

    if limits.shape != (len(properties),2) or len(bfiles) != 2:
        raise ValueError("ndhist limits should have shape (n_properties,2), need 2 bigfiles")

    dset1 = BigFileCatalog(f'{bfiles[0]}',dataset=ptype)
    comm = dset1.comm
    dset2 = BigFileCatalog(f'{bfiles[1]}',dataset=ptype)
    
    #These limits are already in logspace (see below difference)
    edges = np.array([np.linspace(l[0],l[1],num=nbins) for l in limits])
    logger.info(f"edges {edges}")

    #NOTE: I assume here that the snapshots (compressed/not) are in the same order
    #if we want to do some time/model comparison I should change this
    #instead I *can* do a check here to make sure the particles I do sum are the same
    #Of course In the compressed snapshots the particles are in the same order so
    #after checking once I commented this part out
    #pid1 = dset1['ID']
    #pid2 = dset2['ID']
    #same_part = pid1 == pid2
    
    data1 = [parse_property(dset1,prop) for prop in properties]
    data2 = [parse_property(dset2,prop) for prop in properties]
    #a nan property (same_part==0) will not be histogrammed, effectively ignoring it
    #datad = [da.log10(d1/d2 * same_part) for d1,d2 in zip(data1,data2)]
    datad = [da.log10(d1/d2) for d1,d2 in zip(data1,data2)]

    buf,edges = da.histogramdd(datad,edges)

    #since each rank has its own histogram here we can print the total
    #NOTE: a sum less than the list size can be due to a narrow range (limits) OR non-identical IDs

    hist = buf.compute()
    #logger.info(f"Rank {comm.rank} particle list sizes ({pid1.size},{pid2.size}), parts histogrammed {hist.sum()} parts same {same_part.sum().compute()}")

    #each rank has a subset of particles, so summing the histograms will give the final result
    hist = comm.allreduce(hist,op=MPI.SUM)

    return hist,edges
    
def get_diff_1dhist(bfiles,properties,limits,ptype='0/',nbins=100):
    logger.info(f'making log difference 1dhist for {properties}')
    try:
        limits = np.asarray(limits)
    except:
        raise ValueError("ndhist limits should be array-like")

    if limits.shape != (len(properties),2) or len(bfiles) != 2:
        raise ValueError("ndhist limits should have shape (n_properties,2), need 2 bigfiles")

    dset1 = BigFileCatalog(f'{bfiles[0]}',dataset=ptype)
    comm = dset1.comm
    npart = dset1.csize
    dset2 = BigFileCatalog(f'{bfiles[1]}',dataset=ptype)
    
    #These limits are already in logspace (see below difference)
    edges_in = np.array([np.linspace(l[0],l[1],num=nbins) for l in limits])

    #NOTE: see ndhist function comments for future model comparison / ID matching

    hist_arr = np.zeros((len(properties),nbins-1))
    edge_arr = np.zeros((len(properties),nbins))
    for i,prop in enumerate(properties):
        logger.info(f'starting {prop}')
        data1 = parse_property(dset1,prop)
        data2 = parse_property(dset2,prop)
        
        datad = 2*(data2 - data1)/(data2 + data1)
        datad = da.nan_to_num(datad)        
        
        buf,edges = da.histogram(datad,edges_in[i])
        hist = buf.compute()
        
        #each rank has a subset of particles, so summing the histograms will give the final result
        hist = comm.allreduce(hist,op=MPI.SUM)
        hist_arr[i,:] = hist
        edge_arr[i,:] = edges
        logger.info(f"histogrammed {hist.sum()} particles of {npart}")

    return hist_arr,edge_arr


#plot the marginalised histograms from an ND histogram
#CALL ON ONE RANK
def plot_hists(hist,edges,titles,log,nd=False):
    #try to make a nice aspect plot
    comm = MPI.COMM_WORLD
    if comm.rank != 0:
        raise ValueError(f"don't plot on multiple ranks ({comm.rank} of {comm.size})")

    aspect_target = 4./3.
    aspect = 1.
    nrow = 1
    ncol = 1
    hist_iter = range(len(hist.shape) if nd else hist.shape[0])
    for i in hist_iter:
        if i+1 > nrow*ncol:
            if aspect > aspect_target:
                nrow += 1
            else:
                ncol += 1
            aspect = ncol / nrow

    fw = 12
    fh = fw/aspect

    fig,axs = plt.subplots(figsize=(fw, fh),nrows=nrow,ncols=ncol)
    axs = axs.flatten()

    hist_iter = range(len(hist.shape) if nd else len(hist))

    for i in hist_iter:
        #sum over all other axes if n_dimensional
        if nd:
            ax_list = [j for j in range(len(hist.shape))]
            ax_list.remove(i)
            # logger.info(f'axes {i} sum axes {ax_list}')
            # logger.info(f'edges {edges}')
            m_hist = hist.sum(axis=tuple(ax_list))
        else:
            m_hist = hist[i]

        centres = edges[i][:-1] + np.diff(edges[i])/2 if not log[i] else edges[i][:-1] * np.exp(np.diff(np.log(edges[i])))

        if log[i]:
            axs[i].semilogx(centres,m_hist,'k-',linewidth=2)
        else:
            axs[i].plot(centres,m_hist,'k-',linewidth=2)
        axs[i].grid()
        axs[i].set_ylabel('#')
        axs[i].set_xlabel(titles[i])

    return fig,axs

