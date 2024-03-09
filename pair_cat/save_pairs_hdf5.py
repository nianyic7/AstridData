import numpy as np
from bigfile import BigFile
import glob, os, struct

# from astropy.cosmology import FlatLambdaCDM
# import astropy.units as u
# from bh_tools import *
import pickle
import warnings
from scipy.interpolate import interp1d
import argparse


hh = 0.6774
c_mks = 3e8
msun_mks = 2e30
s_to_year = 3.17098e-8
year_to_s = 1.0 / s_to_year
lsun_ergs = 3.9e33
mdot_msun_yr = 1e10 / 980 / 1e6


def calc_lx(mdot):
    """
    input: mdot in Msun/yr
    output: Lx in ergs
    """
    lbol = 0.1 * mdot * msun_mks / year_to_s * c_mks**2
    lbol_lsun = lbol / 3.9e26
    k = 10.83 * (lbol_lsun / 1e10) ** 0.28 + 6.08 * (lbol_lsun / 1e10) ** (-0.02)
    return lbol / k * 1e7


def calc_lbol(mdot):
    """
    input: mdot in Msun/yr
    output: Lx in ergs
    """
    lbol = 0.1 * mdot * msun_mks / year_to_s * c_mks**2
    lbol_ergs = lbol * 1e7
    return lbol_ergs


def edd_ratio(mass, lum):
    return lum / (1.26e38 * mass)


def load_singles(snap, Mmin, Lmin):
    """_summary_

    Args:
        snap (int): snapshot number
        Mmin (float): minimum BH mass in Msun
        Lmin (float): minimum AGN luminosity in erg/s

    Returns:
        np.ndarray: mass, lbol and position of BHs satisfying the above criterion
    """
    print("Loading PIG file: %03d" % snap, flush=True)

    if snap <= 294:
        pig = BigFile("/hildafs/datasets/Asterix/PIG_files/PIG_%03d" % snap)
    else:
        pig = BigFile("/hildafs/datasets/Asterix/PIG2/PIG_%03d" % snap)

    battr = pig["Header"].attrs
    scale_fac = battr["Time"][0]
    redshift = 1.0 / battr["Time"][0] - 1
    print("z=", redshift)

    bhmass = pig.open("5/BlackholeMass")[:] * 1e10 / hh
    bhmdot = pig.open("5/BlackholeAccretionRate")[:].astype(np.float64)
    bhlbol = calc_lbol(bhmdot * mdot_msun_yr)
    bhpos = pig.open("5/Position")[:]
    bhid = pig["5/ID"][:]

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

    f2 = ["id", "mass", "lbol", "pos"]
    dtype2 = ["q", "d", "d", "3d"]
    dt2 = {"names": f2, "formats": dtype2}

    singles = np.zeros(len(mass5), dtype=dt2)
    singles["id"] = id5
    singles["mass"] = mass5
    singles["lbol"] = lbol5
    singles["pos"] = pos5

    return redshift, singles


def pairwise_distance(pos1, pos2):
    dpos = pos1 - pos2
    dpos[dpos > box / 2] -= box
    dpos[dpos < -box / 2] += box
    dr = np.linalg.norm(dpos, axis=-1)
    dr2d = np.sqrt(dpos[:, 0] ** 2 + dpos[:, 1] ** 2)

    return dr, dr2d


def find_duals(singles, rmax, project):
    """_summary_

    Args:
        singles (_type_): _description_
        rmax (_type_): _description_
        project (_type_): _description_

    Returns:
        _type_: _description_
    """
    id5, mass5, lbol5, pos5 = (
        singles["id"],
        singles["mass"],
        singles["lbol"],
        singles["pos"],
    )
    dual_list = []
    Ntot = len(mass5)
    for i in range(Ntot - 1):
        dr, dr2d = pairwise_distance(pos5[i], pos5[i + 1 :, :])
        dr = dr / (1 + redshift) / hh
        dr2d = dr2d / (1 + redshift) / hh
        if project:
            mask = (dr2d <= rmax)
            mask &= (dr2d >= rmin)
            mask = mask.nonzero()[0]
        else:
            mask = (dr <= rmax)
            mask &= (dr >= rmin)
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
                    dr[j-i-1],
                    dr2d[j-i-1],
                )
                for j in js
            ]
    f2 = ["id1", "id2", "m1", "m2", "lb1", "lb2", "pos1", "pos2", "dr", "dr2d"]
    dtype2 = ["q", "q", "d", "d", "d", "d", "3d", "3d", "d", "d"]
    dt2 = {"names": f2, "formats": dtype2}

    duals = np.array(dual_list, dtype=dt2)
    return duals


def find_offsets(singles, rmax, project):
    """_summary_

    Args:
        singles (_type_): _description_
        rmax (_type_): _description_
        project (_type_): _description_

    Returns:
        _type_: _description_
    """
    id5, mass5, lbol5, pos5 = (
        singles["id"],
        singles["mass"],
        singles["lbol"],
        singles["pos"],
    )
    pair_list = []
    Ntot = len(mass5)
    for i in range(Ntot - 1):
        dr, dr2d = pairwise_distance(pos5[i], pos5[i + 1 :, :])
        dr = dr / (1 + redshift) / hh
        dr2d = dr2d / (1 + redshift) / hh
        if project:
            mask = (dr2d <= rmax)
            mask &= (dr2d >= rmin)
            mask = mask.nonzero()[0]
        else:
            mask = (dr <= rmax)
            mask &= (dr >= rmin)
            mask = mask.nonzero()[0]

        if len(mask) >= 1:
            js = mask + i + 1  # abs index
            pair_list += [
                (
                    id5[i],
                    id5[j],
                    mass5[i],
                    mass5[j],
                    lbol5[i],
                    lbol5[j],
                    pos5[i],
                    pos5[j],
                    dr[j-i-1],
                    dr2d[j-i-1],
                )
                for j in js
            ]
    f2 = ["id1", "id2", "m1", "m2", "lb1", "lb2", "pos1", "pos2", "dr", "dr2d"]
    dtype2 = ["q", "q", "d", "d", "d", "d", "3d", "3d", "d", "d"]
    dt2 = {"names": f2, "formats": dtype2}

    pairs = np.array(pair_list, dtype=dt2)
    # lum cut in the end
    print("total pairs before lum cut:", len(pairs), flush=True)
    L_bright = np.maximum(pairs["lb1"], pairs["lb2"])
    L_faint = np.minimum(pairs["lb1"], pairs["lb2"])
    mask = L_bright >= Lmin
    mask &= L_faint < Lmin
    offsets = pairs[mask]
    
    return offsets


if __name__ == "__main__":
    # ------------------- Arguments ------------------------------------
    parser = argparse.ArgumentParser(description="save-duals")
    parser.add_argument("--snap", required=True, type=int, help="index of the PIG file")
    parser.add_argument(
        "--rmax", required=True, type=float, help="maximum separation for pair in kpc"
    )
    parser.add_argument(
        "--rmin", required=True, type=float, help="minimum separation for pair in kpc"
    )
    parser.add_argument(
        "--Mmin",
        default=0,
        type=float,
        help="minimum mass in Msun, no cut if not specified",
    )
    parser.add_argument(
        "--Lmin",
        default=0,
        type=float,
        help="minimum Lbol in erg/s, no cut if not specified",
    )
    parser.add_argument(
        "--project",
        action="store_true",
        help="whether to use projected distance, default to False if not flagged",
    )
    parser.add_argument(
        "--offset",
        action="store_true",
        help="whether to find offset, default to False if not flagged",
    )
    parser.add_argument(
        "--savesingle",
        action="store_true",
        help="save singles satisfying the selection criterion, default to False if not flagged",
    )

    args = parser.parse_args()
    snap = int(args.snap)
    rmax = args.rmax
    rmin = args.rmin
    Mmin = args.Mmin
    Lmin = args.Lmin
    offset = args.offset
    project = args.project
    save_single = args.savesingle
    # ------------------ Read in BHs -----------------------------------
    box = 250000.0
    if offset:
        # no lum cut for singles
        redshift, singles = load_singles(snap=snap, Mmin=Mmin, Lmin=0)
        pairs = find_offsets(singles, rmax, project)
        print("Number of offsets:", len(pairs), flush=True)
    else:
        redshift, singles = load_singles(snap=snap, Mmin=Mmin, Lmin=Lmin)
        pairs = find_duals(singles, rmax, project)
        print("Number of duals:", len(pairs), flush=True)

        
    

    # swap sp that m1 1 is the massive one
    swap = pairs["m1"] < pairs["m2"]
    pairs["m1"][swap], pairs["m2"][swap] = pairs["m2"][swap], pairs["m1"][swap]
    pairs["id1"][swap], pairs["id2"][swap] = pairs["id2"][swap], pairs["id1"][swap]
    pairs["lb1"][swap], pairs["lb2"][swap] = pairs["lb2"][swap], pairs["lb1"][swap]
    if offset:
        path = "/hildafs/home/nianyic/scratch1/dualAGN/offset_catalog_yueshen/"
    else:
        path = "/hildafs/home/nianyic/scratch1/dualAGN/dual_catalog_yueshen/"
    
    if project:
        np.save(
            path + "pair-s%d-m%.1e-lb%.1e-dr2d=%d-%dkpc.npy" % (snap, Mmin, Lmin, rmin, rmax),
            pairs,
        )
    else:
        np.save(
            path + "pair-s%d-m%.1e-lb%.1e-dr3d=%d-%dkpc.npy" % (snap, Mmin, Lmin, rmin, rmax),
            pairs,
        )
    if save_single:
        np.save(path + "singles-s%d-m%.1e-lb%.1e.npy" % (snap, Mmin, Lmin), singles)

    print("Done!", flush=True)
