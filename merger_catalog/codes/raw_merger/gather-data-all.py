import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from bigfile import BigFile
import glob,os,struct
from scipy.interpolate import interp1d
from colossus.cosmology import cosmology
import seaborn as sns

sns.set()
sns.set_palette("Set2")
sns.set_style('ticks',{'ytick.direction':'in','xtick.direction':'in'})

cmap = plt.get_cmap("Set2")
sns.set_context("paper", font_scale=1.7,rc={"axes.linewidth": 1.3,"lines.linewidth": 2.5,"patch.linewidth": 2.2})
from matplotlib import rcParams as rc
import matplotlib as mpl

def load_astrid_cosmo():
    pig = BigFile('/hildafs/datasets/Asterix/PIG_files/PIG_214')
    battr = pig["Header"].attrs
    Lbox = battr['BoxSize']
    hh = battr['HubbleParam']
    om0 = battr['Omega0']
    omb = battr['OmegaBaryon']
    oml = battr['OmegaLambda']
    Nfof = battr['NumFOFGroupsTotal']
    sigma8 = 0.82

    # set-up cosmology
    params = {'flat': True, 'H0': 100*hh, 'Om0': om0, 'Ob0': omb, 'sigma8': sigma8, 'ns': 0.96}
    cosmo = cosmology.setCosmology('myCosmo', params)
    MassTable = battr["MassTable"]
    
    return Lbox, MassTable, cosmo

snaps = np.loadtxt('/hildafs/datasets/Asterix/PIG2/Snapshots.txt')
snap2z = {s[0]:1./s[1]-1 for s in snaps}


# some constants and unit conversions
msun_mks = 1.989e30
pc_mks = 3.086e16
grav_mks = 6.67e-11
km_mks = 1e3
yr_mks = 3.154e+7
c_mks = 3e8

Mpc_to_m = 3.086e+22
m_to_Mpc = 1./Mpc_to_m
s_to_year = 3.17098e-8
c_Mpc_yr = c_mks*m_to_Mpc/s_to_year
fac = 121.14740013761634  # GADGET unit Protonmass / Bolztman const
GAMMA = 5 / 3.
Xh=0.76
hh = h = 0.6774

Lbox, MassTable, cosmo = load_astrid_cosmo()
# conversion between time and redshift
z_arr = np.linspace(0,20,1000)
time = cosmo.age(z_arr) # Gyr
def z_to_t(x):
    return interp1d(z_arr, time,fill_value='extrapolate')(x)
def t_to_z(x):
    return interp1d(time, z_arr,fill_value='extrapolate')(x) 


merger_root = "/hildafs/datasets/Asterix/BH-mergers"
merger_files = glob.glob(merger_root + "/bh-merger-extended*.npy")
merger_files.sort()

mergers = []
for file in merger_files:
    mergers.append(np.load(file))
mergers = np.concatenate(mergers)

print("merger data fields:", mergers.dtype.names)
print("merger shape:", mergers.shape)

swbh = np.minimum(mergers["ID1"], mergers["ID2"])

uid, uidx = np.unique(swbh, return_index=True)
mergers = mergers[uidx]

print("unique mergers:", len(mergers))
savedir = "/hildafs/home/nianyic/scratch1/Astrid_data/lisa_mbhcat/up_to_date_mergers_nogal.npy"
np.save(savedir, mergers)
print("saved catalog to:", savedir)


fig, ax = plt.subplots(1, 1, figsize=(6, 5))
vsim = 250**3
zmin = min(mergers["z"])
zmax = max(mergers["z"])
print(f"zmin={zmin}, zmax={zmax}")
zbin = np.linspace(zmin,zmax,20)
dz = zbin[1:] - zbin[:-1]
zc = 0.5*(zbin[1:]+zbin[:-1])
rc = np.array([cosmo.comovingDistance(0.,zz) for zz in zc])

n0,e0 = np.histogram(np.concatenate([mergers['z']]),bins=zbin)

rate0 = n0/dz/vsim*4*np.pi*rc**2*c_Mpc_yr

ax.plot(zc, rate0, label="all")
ax.grid()
ax.set(xlabel="redshift", ylabel=r"dN/dzdt$\rm{[Myr]^{-1}}$")
plt.savefig("merger_rate_uptodate.png", bbox_inches="tight", dpi=130)
print("saved merger rate plot")
