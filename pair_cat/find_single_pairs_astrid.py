import numpy as np
from bigfile import BigFile
import glob, os, struct
import warnings
import datetime
from utils import *


hh = 0.6774


IntTypes = set(
    [
        np.dtype("uint64"),
        np.dtype("uint32"),
        np.dtype("uint8"),
        np.dtype("int64"),
        np.dtype("int32"),
        np.dtype("int8"),
    ]
)


KeepCols = set(["ElectronAbundance", "Metallicity", "Metals", "BlackholeSwallowTime"])

def recover(col, odata):
    """
    Recover compressed data into original datatypes and values

    Args:
        col (str): data column name, e.g. "Position"
        odata (np.ndarray): array to be recovered

    Returns:
        np.ndarray: array recovered from compressed data
    """
    if odata.dtype in IntTypes:
        return odata
    elif col == "Position":
        rdata = odata.astype(np.dtype("float64"))
        return rdata
    elif col in KeepCols:
        return odata
    else:
        rdata = odata.astype(np.dtype("float32"))
        mask = rdata > 0
        rdata[mask] = np.power(rdata[mask], 10)
        mask = rdata < 0
        rdata[mask] = -np.power(rdata[mask], 10)

    return rdata



def set_metadata_astrid():
        metadata = {
        # ---- Header data - identifying information for the dataset
        # REQUIRED
        'SimulationName': 'Astrid',
        'SimulationVersion': '1.0',
        'ModelType': 'Hydro',
        'Version': '1.0',
        'Date': str(datetime.datetime.now()),
        'Contributor': ["Nianyi Chen", "Yueying Ni"],
        'Email': ["nianyic@andrew.cmu.edu", "yueying.ni@cfa.harvard.edu"],
        'Principal': ["Tiziana Di Matteo"],
        'Reference': ["10.1093/mnras/stac1432", "10.1093/mnras/stac351", "10.1093/mnras/stad3084"],
        # OPTIONAL:
        'Website': ["https://astrid-portal.psc.edu"],

        # ---- Model parameters - metadata specification for simulation(s) used to construct catalog
        # REQUIRED
        'HubbleConstant': 67.74, #km s^-1 Mpc^-1
        'BoxSize': 369, # cMpc
        'MinBHSeedMass': 4.4e4, #M_sun
        'MinRedshift': 1.316,
        'MaxRedshift': 99,
        'StellarMassResolution': 4.7e5, # M_sun
        'DarkMatterMassResolution': 9.63e6, # M_sun
        'SpatialResolution': np.nan, # kpc
        # OPTIONAL:
        'MinimumDarkMatterHaloMass': 3e8, # M_sun
        'GasMassResolution': 1.87e6, # M_sun
        }
        return metadata


def load_singles_astrid(datadir, snap, Mmin, Lmin):
    """_summary_

    Args:
        snap (int): snapshot number
        Mmin (float): minimum BH mass in Msun
        Lmin (float): minimum AGN luminosity in erg/s

    Returns:
        dictionary:
            redshift (float): redshift of the snapshot
            snap (int): snapshot number
            boxsize (float): boxsize in kpc/h
            bhmass (np.array): BH mass in Msun
            bhmdot (np.array): BH accretion rate in Msun/yr
            bhlbol (np.array): BH bolometric luminosity in erg/s
            bhpos (np.array): BH position in kpc/h
            bhid (np.array): BH ID
    """
    print("Loading PIG file: %03d" % snap, flush=True)


    pig = BigFile(datadir + "/PIG_%03d" % snap)
    battr = pig["Header"].attrs
    scale_fac = battr["Time"][0]
    boxsize = battr["BoxSize"][0]
    redshift = 1.0 / battr["Time"][0] - 1
    hubble = hh = battr["HubbleParam"][0]
    print("z=", redshift)

    # check if file is compressed:
    if pig["5/BlackholeMass"].dtype == np.dtype("float16"):
        print("Compressed file detected, recovering data...", flush=True)
        bhmass = recover("BlackholeMass", pig["5/BlackholeMass"][:])
        bhmass *= 1e10 / hh
        bhmdot = recover("BlackholeAccretionRate", pig["5/BlackholeAccretionRate"][:])
        bhmdot = bhmdot.astype(np.float64)
        bhpos = recover("Position", pig["5/Position"][:])
    else:
        bhmass = pig["5/BlackholeMass"][:] * 1e10 / hh
        bhmdot = pig["5/BlackholeAccretionRate"][:].astype(np.float64)
        bhpos = pig["5/Position"][:]

    bhid = pig["5/ID"][:]
    bhgid = pig["5/GroupID"][:]
    bhlbol = calc_lbol(bhmdot * mdot_msun_yr)

    print(
        "Selecting BH with M > %.1e Msun, Lbol > %.1e erg/s" % (Mmin, Lmin), flush=True
    )

    maskAGN = bhmass > Mmin
    maskAGN &= bhlbol > Lmin

    bidx = maskAGN.nonzero()[0]
    print("Seleciton resulting in %d AGNs" % len(bidx), flush=True)

    mass5 = bhmass[bidx]
    lbol5 = bhlbol[bidx]
    pos5 = bhpos[bidx]
    id5 = bhid[bidx]
    gid5 = bhgid[bidx]
    mdot = bhmdot[bidx] * mdot_msun_yr

    f2 = ["id", "mass", "lbol", "pos"]
    dtype2 = ["q", "d", "d", "3d"]
    dt2 = {"names": f2, "formats": dtype2}

    singles = {}
    singles["bhmdot"] = mdot
    singles["bhid"] = id5
    singles["bhmass"] = mass5
    singles["bhlbol"] = lbol5
    singles["bhpos"] = pos5
    singles["bhgid"] = gid5
    singles["redshift"] = redshift
    singles["snap"] = snap
    singles["boxsize"] = boxsize
    singles["hubble"] = hubble
    return singles


def pairwise_distance(pos1, pos2, box):
    dpos = pos1 - pos2
    dpos[dpos > box / 2] -= box
    dpos[dpos < -box / 2] += box
    dr = np.linalg.norm(dpos, axis=-1)
    dr2d = np.sqrt(dpos[:, 0] ** 2 + dpos[:, 1] ** 2)

    return dr, dr2d


def find_pairs_astrid(singles, rmax):
    """_summary_

    Args:
        singles (dict): single bh catalog returned by load_singles
        rmax (float): maximum separation in kpc for pair search

    Returns:
        dictionary:
            redshift (float): redshift of the snapshot
            snap (int): snapshot number
            bhmass (np.array): BH mass in Msun
            bhmdot (np.array): BH accretion rate in Msun/yr
            bhlbol (np.array): BH bolometric luminosity in erg/s
            bhpos (np.array): BH position in kpc/h
            bhid (np.array): BH ID

    """
    id5, mass5, lbol5, pos5 = (
        singles["bhid"],
        singles["bhmass"],
        singles["bhlbol"],
        singles["bhpos"],
    )
    gid5, mdot5 = singles["bhgid"], singles["bhmdot"]
    redshift = singles["redshift"]
    boxsize = singles["boxsize"]

    dual_list = []
    Ntot = len(mass5)
    # Ntot = 5000
    for i in range(Ntot - 1):
        dr, dr2d = pairwise_distance(pos5[i], pos5[i + 1 :, :], boxsize)
        dr = dr / (1 + redshift) / hh
        dr2d = dr2d / (1 + redshift) / hh

        mask = (dr <= rmax)
        mask = mask.nonzero()[0]
        if len(mask) >= 1:
            js = mask + i + 1  # abs index
            dual_list += [
                (
                    id5[i],
                    id5[j],
                    mass5[i],
                    mass5[j],
                    lbol5[i],
                    lbol5[j],
                    pos5[i],
                    pos5[j],
                    gid5[i],
                    gid5[j],
                    mdot5[i],
                    mdot5[j],
                    dr[j-i-1],
                    dr2d[j-i-1],
                )
                for j in js
            ]

    
    f2 = ["id1", "id2", "m1", "m2", "lb1", "lb2", "pos1",\
         "pos2", "gid1", "gid2", "mdot1", "mdot2", "dr", "dr2d"]
    dtype2 = ["q", "q", "d", "d", "d", "d", "3d", "3d", \
        "q", "q", "d", "d", "d", "d"]
    dt2 = {"names": f2, "formats": dtype2}
    pairs = np.array(dual_list, dtype=dt2)

    pairs_dict = {}
    for k in f2:
        pairs_dict[k] = pairs[k]
    pairs_dict["redshift"] = redshift
    pairs_dict["snap"] = singles["snap"]
    pairs_dict["boxsize"] = boxsize
    pairs_dict["hubble"] = singles["hubble"]
    pairs_dict["rmax"] = rmax
    print("found %d pairs"%len(pairs))
    return pairs_dict