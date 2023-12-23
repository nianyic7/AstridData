"""
extract merger event from ASTRID bhdetail files in bigfile format
"""
import numpy as np
from bigfile import BigFile
import struct
import glob
import os, sys
import argparse

hh = 0.6774


def id2_by_swallowed(zmerge, id1):
    sidx = (SWID == id1).nonzero()[0]
    # find the id2 swallowed closest to zmerge
    id2 = None
    zdiff_min = 0.01
    for ii in ID2[sidx]:
        mask2 = BHID == ii
        zorder = np.argsort(-redshift[mask2])
        z2 = redshift[mask2][zorder]
        sw2 = Swallowed[mask2][zorder]
        swidx = (sw2 > 0).nonzero()[0][0]
        zsw = z2[swidx]
        if np.abs(zsw - zmerge) < zdiff_min:
            zdiff_min = np.abs(zsw - zmerge)
            id2 = ii
    print(
        "found ID2=%d swallowed at %.4f from from zmerge=%.4f"(id2, zdiff_min, zmerge),
        flush=True,
    )
    return id2


def id2_by_mass_separation(zmerge, id1, pos1, tarmass, zwindow=1e-2):
    zmask = np.abs(redshift - zmerge) < zwindow
    idx = np.where(BHMass[zmask] == tarmass)[0]
    id2 = None

    if len(idx) > 0:
        for i in idx:
            pos2 = BHpos[zmask][i]
            rr = np.linalg.norm(pos2 - pos1)
            if rr < 5:
                id2 = BHID[zmask][i]
                print("found ID2 = ", id2, "r = %.2f" % rr)
    else:
        print(
            "within dz=%.5f,didn't find matched ID2 with acBHmass" % zwindow, flush=True
        )
    return id2


def id2_by_swid(zmerge, id1):
    sidx = (SWID == id1).nonzero()[0]
    # -----------------------------------------
    if len(sidx) == 0:
        return None
    elif len(sidx) == 1:
        id2 = ID2[sidx[0]]  # a unique swallowee
        return id2
    else:
        # id1 swallowed more than one bhs
        # need to match redshift
        for id2 in ID2[sidx]:
            mask2 = BHID == id2
            z2 = redshift[mask2]
            sw2 = swallowed[mask2]
            zmask = z2 <= zmerge
            zmask &= z2 >= zend
            if (np.all(sw2[z2 > zmerge] == 0)) and (np.all(sw2[zmask] == 1)):
                return id2
    return None


def find_merger(zwindow=1e-3):
    # search by swallower
    mask = acBHMass > 1e-8
    mask &= Swallowed == 0
    acidx = mask.nonzero()[0]
    print("expected mergers from acBHMass", len(acidx))
    if len(swidx) == 0:
        print("No merger in this chunk", flush=True)
        return

    for idx in swidx:
        zmerge = redshift[idx]
        id1 = BHID[idx]  # swallower
        id2 = id2_by_swid(zmerge, id1)
        if id2 is None:
            # second pass
            print(
                "zmerge=%.3f, id1=%d, acBHMass=%.7f, swid matching failed, searching by mass and separation..."
                % (zmerge, id1, acBHMass[idx]),
                flush=True,
            )
            id2 = id2_by_mass_separation(
                zmerge, id1, pos1=BHpos[idx], tarmass=acBHMass[idx]
            )
        if id2 is None:
            print(
                "zmerge=%.3f, id1=%d, acBHMass=%.7f, position matching failed, searching by swallow redshift..."
                % (zmerge, id1, acBHMass[idx]),
                flush=True,
            )
            # third pass
            id2 = id2_by_swallowed(zmerge, id1)
        if id2 is None:
            print("still did not find matched a BH2!!!", "acBHMass = ", acBHMass[idx])
            continue
        # find masses
        m1 = BHMass[idx] - acBHMass[idx]
        m2 = acBHMass[idx]
        if m1 < m2:  # convention is m1>m2
            m1, m2 = m2, m1
            id1, id2 = id2, id1

        mergers.append((zmerge, id1, id2, m1 * 1e10 / hh, m2 * 1e10 / hh))

    f1 = ["z", "ID1", "ID2", "m1", "m2"]
    d1 = ["d", "q", "q", "d", "d"]
    dt = {"names": f1, "formats": d1}
    mergers = np.array(mergers, dtype=dt)
    return mergers


mdata = get_merger(bh_reduce_file)
if mdata is not None:
    mdata = np.unique(mdata)
    print("unique event number:", len(mdata))
    np.save(ofilename, mdata)


if __name__ == "__main__":
    # -------------- Cmd line Args ------------------------------
    parser = argparse.ArgumentParser(description="find mergers")
    parser.add_argument(
        "--indir",
        required=True,
        type=str,
        help="path of the bh detail bigfile",
    )
    parser.add_argument(
        "--snap", required=True, type=int, help="bigfile snapshot number"
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=str,
        help="path of the output merger data directory",
    )
    # ------------------------- args -----------------------------------
    args = parser.parse_args()
    snap = str(args.snap)
    bhfile = indir + "/BH-Details-R" + snap
    ofile = outdir + "/mdata/bh-merger-R" + snap
    print("reading", bhfile, flush=True)

    # ------------------------- load data -----------------------------------
    bf = BigFile(bhfile)
    acBHMass = bf["acBHMass"][:]
    Swallowed = bf["Swallowed"][:]
    SwallowID = bf["SwallowID"][:]
    BHID = bf["BHID"][:]
    redshift = bf["z"][:]
    BHpos = bf["BHpos"][:]
    BHMass = bf["BHMass"][:]
    zstart, zend = np.max(redshift), np.min(redshift[redshift > 0])
    print("redshift range of this chunk:", zstart, zend, flush=True)

    # find SWID
    mask = SwallowID > 0
    SWID = SwallowID[mask]  # swallower
    ID2 = BHID[mask]  # swallowee

    ID2, idx = np.unique(ID2, return_index=True)
    SWID = SWID[idx]

    print("Unique SwallowID-BHID pairs", len(SWID))
    print("Unique SwallowID", len(np.unique(SWID)), flush=True)

    # ------------------------ find mergers ---------------------------
    mergers = find_merger(zwindow=1e-3)

    if mergers is not None:
        np.save(ofile, mergers)
    print("Saved data to:", ofile, flush=True)
