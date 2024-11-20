"""
extract merger event from bhdetail-reduce files
Usage: python extract-merger-v1.py $snap

An extension of extract-merger-v1.py,
contains additional information about mergers:
positions and accretion rates just before the merger.
"""

import numpy as np
from bigfile import BigFile
import sys
import glob


hh = 0.6774
mdot_msun_yr = 1e10/978/1e6
snap = None
cluster = None

acBHMass_ds = None
swallowed_ds = None
swallowID_ds = None
BHID_ds = None
redshift_ds = None
BHpos_ds = None
BHvel_ds = None
BHMdot_ds = None
BHMass_ds = None
zstart = None
zend = None

# prepare for the previous chunk
BHID_pr = None
redshift_pr = None

BHpos_pr = None
BHvel_pr = None
BHMass_pr = None
BHMdot_pr = None

cur_chunk_i = 0



def load_bh_data(filename):
    global acBHMass_ds, swallowed_ds, swallowID_ds, BHID_ds, redshift_ds, \
        BHpos_ds, BHvel_ds, BHMass_ds, BHMdot_ds, zstart, zend

    print(f"Reading {filename}")
    bhd = BigFile(filename)
    # load partial data to check if there is merger event
    acBHMass_ds = bhd.open('acBHMass')[:]
    swallowed_ds = bhd.open('Swallowed')[:].astype(np.bool)
    # check if no merger event
    mask0 = (acBHMass_ds > 1e-8) & (swallowed_ds == 0)
    if np.sum(mask0) < 1:
        sys.exit("No merger event found.")

    # load full data
    BHID_ds = bhd.open('BHID')[:]
    redshift_ds = bhd.open('z')[:]
    # sort by ID then redshift (high to low)
    sort_indices = np.lexsort((-redshift_ds, BHID_ds))
    BHID_ds = BHID_ds[sort_indices]
    redshift_ds = redshift_ds[sort_indices]
    acBHMass_ds = acBHMass_ds[sort_indices]
    swallowed_ds = swallowed_ds[sort_indices]

    swallowID_ds = bhd.open('SwallowID')[:][sort_indices]
    BHpos_ds = bhd.open('BHpos')[:].astype(np.float32)
    BHpos_ds = BHpos_ds[sort_indices]
    BHvel_ds = bhd.open('BHvel')[:][sort_indices]
    BHMass_ds = bhd.open('BHMass')[:][sort_indices]
    BHMdot_ds = bhd.open('Mdot')[:][sort_indices]

    # get redshift range
    zstart, zend = np.max(redshift_ds), np.min(redshift_ds[redshift_ds > 0])
    print("load_bh_data, Total length:", len(BHID_ds))
    print("Redshift range:", zstart, zend, flush=True)



def get_previous_chunk(snap):
    global BHID_pr, redshift_pr, \
        BHpos_pr, BHvel_pr, BHMass_pr, BHMdot_pr

    fidx_to_path, snap_to_fidx = get_all_details_file()
    curr_idx = snap_to_fidx[snap]
    prev_path = fidx_to_path[curr_idx - 1]
    bf_prev = BigFile(prev_path)

    print(f"Reading previous chunk {prev_path}", flush=True)

    # load data
    BHID_pr = bf_prev.open('BHID')[:]
    redshift_pr = bf_prev.open('z')[:]
    zmin = np.min(redshift_pr[redshift_pr > 0])
    print(f"Redshift range of the previous chunk: {np.max(redshift_pr)}, {zmin}", flush=True)
    mask = redshift_pr < zmin + 0.02
    BHID_pr = BHID_pr[mask]
    redshift_pr = redshift_pr[mask]
    # sort by ID then redshift (high to low)
    sort_indices = np.lexsort((-redshift_pr, BHID_pr))
    BHID_pr = BHID_pr[sort_indices]
    redshift_pr = redshift_pr[sort_indices]

    BHpos_pr = bf_prev.open('BHpos')[:][mask].astype(np.float32)
    BHpos_pr = BHpos_pr[sort_indices]
    BHvel_pr = bf_prev.open('BHvel')[:][mask][sort_indices]
    BHMass_pr = bf_prev.open('BHMass')[:][mask][sort_indices]
    BHMdot_pr = bf_prev.open('Mdot')[:][mask][sort_indices]



def get_acBH_events():
    """_summary_

    Returns:
        _type_: _description_
    """
    # These are indices immediately after the mergers
    mask0 = (acBHMass_ds > 1e-8) & (swallowed_ds == 0)
    acBH_indices = mask0.nonzero()[0]
    Id1s = BHID_ds[acBH_indices]
    print(f"Number of mergers from acBHMass: {len(acBH_indices)}")

    index1_to_ID1 = {acBH_indices[i]: Id1 for i, Id1 in enumerate(Id1s)}
    return index1_to_ID1

def get_swallowed_events():
    """_summary_

    Returns:
        _type_: _description_
    """
    # get ids
    mask1 = (swallowID_ds > 0).nonzero()[0]
    ID1 = swallowID_ds[mask1]  # all swallower
    ID2 = BHID_ds[mask1]  # the one got swallowed

    # one BH can not be swallowed twice, so reduce ID2
    # np unique always returns the first occurrence, which is nice
    ID2, i = np.unique(ID2, return_index=True)
    ID1 = ID1[i]
    mask1 = mask1[i]
    
    print("Find unique SwallowID-BHID pairs", len(ID1))
    print("Number of unique SwallowID", len(np.unique(ID1)), flush=True)

    ID2_to_index2 = {Id2: mask1[i] for i, Id2 in enumerate(ID2)}
    ID2_to_ID1 = {Id2: ID1[i] for i, Id2 in enumerate(ID2)}
    return ID2_to_index2, ID2_to_ID1


def get_index_data():
    index1_to_ID1 = get_acBH_events() # indices of BH1 at merger
    # indices of BH2 at merger, note: this may be a subset of ID1_to_index1
    ID2_to_index2, ID2_to_ID1 = get_swallowed_events() 

    idx1_ac = np.array(list(index1_to_ID1.keys()))
    ID1_sw = np.array(list(ID2_to_ID1.values()))
    ID2_sw = np.array(list(ID2_to_index2.keys()))

    num_mergers = 0
    data = []
    case1, case2, case3, case4, case5 = 0, 0, 0, 0, 0


    for idx1 in idx1_ac:
        id_1 = index1_to_ID1[idx1]
        assert id_1 == BHID_ds[idx1]

        zmerge = redshift_ds[idx1]
        # all events where ID1 is swallower
        swallowed_idx = np.where(ID1_sw == id_1)[0]

        id_2 = None
        flag = None
        if len(swallowed_idx) == 0: 
            # case1: BH1 swallowed the mass, but is not recoreded on BH2 history
            print(f"ID1={id_1}, Mass={BHMass_ds[idx1]}, acBHMass={acBHMass_ds[idx1]}, "
                  f"now trying to find ID2 by position and mass", flush=True)
            id_2, idx2 = find_ID2_by_pos_mass(zmerge, acBHMass_ds[idx1], BHpos_ds[idx1])
            if id_2 is None:
                print(f"Still did not find matched a BH2, !!! acBHMass = {acBHMass_ds[idx1]}", flush=True)
                case5 += 1
                continue
            else:
                case1 += 1
                flag = 1
                data.append((id_1, id_2, idx1, idx2, flag))

        elif len(swallowed_idx) == 1:
            # case2: the easy case, only one BH2 is swallowed and correctly recorded
            id_2 = ID2_sw[swallowed_idx[0]]
            idx2 = ID2_to_index2[id_2]
            case2 += 1
            flag = 2
            data.append((id_1, id_2, idx1, idx2, flag))

        else:
            # case3: multiple BH2 are swallowed, find the one logged and closest to z_acBH
            for potential_id in ID2_sw[swallowed_idx]:
                if validate_ID2_by_z(potential_id, zmerge):
                    id_2 = potential_id
                    idx2 = ID2_to_index2[id_2]
                    case3 += 1
                    flag = 3
                    data.append((id_1, id_2, idx1, idx2, flag))
                    # no break because there may be two swallows at zmerge

        # if still no id_2, try to find it by zflip
        if ((len(swallowed_idx) > 1) & (id_2 is None)):
            zflip_idx2 = np.array([find_z_when_swallowed_flips(potential_id)
                            for potential_id in ID2_sw[swallowed_idx]])
            z_flip = zflip_idx2[:, 0]
            idx2 = zflip_idx2[:, 1]

            closest_z_flip_idx = np.argmin(np.abs(z_flip - zmerge))
            id_2 = ID2_sw[swallowed_idx[closest_z_flip_idx]]
            idx2 = idx2[closest_z_flip_idx]
            print(
                f"With second zflip search, find ID2 = {id_2}, zflips = {z_flip}", flush=True)
            case4 += 1
            flag = 4
            data.append((id_1, id_2, idx1, idx2, flag))

        if id_2 is None:
            case5 += 1
            print(
                f"Still did not find matched a BH2, !!! acBHMass = {acBHMass_ds[idx1]}", flush=True)
            # continue
    dt = np.dtype([('ID1', 'q'), ('ID2', 'q'), ('idx1', 'i'), ('idx2', 'i'), ('flag', 'i')])
    print("Total number of mergers found:", len(data), flush=True)
    print(f"case1: {case1}, case2: {case2}, case3: {case3}, case4: {case4}, case5: {case5}", flush=True)

    data = np.array(data, dtype=dt)
    return data

    



def get_bh_info(events):
    """Given all event index, get the corresponding BH1 and BH2 info

    Args:
        events (np.array):
            [('ID1', 'q'), ('ID2', 'q'), ('idx1', 'i'), ('idx2', 'i'), ('flag', 'i')]

    Returns:
        merger data:
        [('z', 'd'), ('ID1', 'q'), ('ID2', 'q'), ('m1', 'd'), ('m2', 'd'), ('mdot1', 'd'), \
            ('mdot2', 'd'), ('pos1', '3d'), ('pos2', '3d'), ('v1', '3d'), ('v2', '3d')])
    """
    # initialize empty np array for merger_data
    merger_dtype = np.dtype(
        [('z', 'd'), ('ID1', 'q'), ('ID2', 'q'), ('m1', 'd'), ('m2', 'd'), ('mdot1', 'd'), \
            ('mdot2', 'd'), ('pos1', '3d'), ('pos2', '3d'), ('v1', '3d'), ('v2', '3d')])
    merger_data = []
    for i, (id1, id2, idx1, idx2, flag) in enumerate(events):
        zmerge = redshift_ds[idx1]
        m1 = BHMass_ds[idx1] - acBHMass_ds[idx1]
        m2 = acBHMass_ds[idx1]
        id1_1 = BHID_ds[idx1 - 1]
        id2_1 = BHID_ds[idx2 - 1]

        #-----------------------------------------------------------------------------------
        if id1_1 != id1 or id2_1 != id2:
            print("swallow happened at the very beginning of the chunk,", flush=True)
            print(f"ID1 = {id1}, ID1_1 = {id1_1}, ID2 = {id2}, ID2_1 = {id2_1}", flush=True)
            print(f"zmerge = {zmerge}, zstart = {zstart}", flush=True)
            if BHID_pr is None:
                try:
                    get_previous_chunk(snap)
                except:
                    print("No previous chunk, discard this event", flush=True)
                    continue

                # get_previous_chunk(snap)

            mask1 = (BHID_pr == id1).nonzero()[0]
            mask2 = (BHID_pr == id2).nonzero()[0]
            if len(mask1) == 0:
                print(f"ID1 = {id1} not found in the previous chunk, discard this event", flush=True)
                continue
            if len(mask2) == 0:
                print(f"ID2 = {id2} not found in the previous chunk, discard this event", flush=True)
                continue

            idx1_1 = mask1[-1]
            idx2_1 = mask2[-1]
            # pointer to the previous chunk
            BHID, redshift, BHpos, BHvel, BHMass, BHMdot \
            = BHID_pr, redshift_pr, BHpos_pr, BHvel_pr, BHMass_pr, BHMdot_pr
        else:
            idx1_1 = idx1 - 1
            idx2_1 = idx2 - 1
            # pointer to the current chunk
            BHID, redshift, BHpos, BHvel, BHMass, BHMdot \
            = BHID_ds, redshift_ds, BHpos_ds, BHvel_ds, BHMass_ds, BHMdot_ds
        #------------------------------------------------------------------------------------  

        while redshift[idx1_1] <= zmerge:
            idx1_1 -= 1
        while redshift[idx2_1] <= zmerge:
            idx2_1 -= 1
        m1_alt = BHMass[idx1_1]
        m2_alt = BHMass[idx2_1]
        z_prev1 = redshift[idx1_1]
        z_prev2 = redshift[idx2_1]
        zmerge1 = redshift_ds[idx1]
        zmerge2 = redshift_ds[idx2]
        if not np.isclose(m1, m1_alt, rtol=3e-2):
            print(f"Swallower mass info from before/after merger differs, record the current timestep m1 \
            m1 = {m1}, m1_alt = {m1_alt}, m2 = {m2}, m2_alt = {m2_alt}, \
                z_prev1 = {z_prev1}, z_prev2 = {z_prev2}, \
                    zmerge1 = {zmerge1}, zmerge2 = {zmerge2}, flag = {flag}")
        if not np.isclose(m2, m2_alt, rtol=3e-2):
            print(f"A case of a simultaneous swallow,\
                let's just record the previous timestep m2", flush=True)
            m2 = m2_alt
        # convention is m1 > m2
        if m1 < m2:
            m1, m2 = m2, m1
            id1, id2 = id2, id1
            idx1, idx2 = idx2, idx1
        # get mdot info
        mdot1 = BHMdot[idx1_1] * mdot_msun_yr # to Msun/yr
        mdot2 = BHMdot[idx2_1] * mdot_msun_yr # to Msun/yr
        pos1 = BHpos[idx1_1] # in ckpc/h
        pos2 = BHpos[idx2_1] # in ckpc/h
        vel1 = BHvel[idx1_1] * (1 + zmerge) # to km/s
        vel2 = BHvel[idx2_1] * (1 + zmerge) # to km/s

        merger_data.append((zmerge, id1, id2, m1 * 1e10 / hh, m2 * 1e10 / hh, \
            mdot1, mdot2, pos1, pos2, vel1, vel2))
            
    #----------------------------------------------------------------------------------------------
    merger_data = np.array(merger_data, dtype=merger_dtype)
    return merger_data


def find_ID2_by_pos_mass(zmerge, target_acBHMass, target_BHpos):
    """
    If the swallowed timestep is not logged, we may find ID2 via BHpos, BHmass
    Use chunk for faster processing. Since dataset is ordered by redshift,
    we can start from the last chunk we processed
    """
    global cur_chunk_i

    # process in chunks to avoid longer time
    chunk_size = 2**17  # 1000000 # 2**17 for RM, 1e6 for HENON?
    num_chunks = len(redshift_ds) // chunk_size + 1

    id2 = None
    idx2 = None
    for i in range(cur_chunk_i, num_chunks):
        # get start/end index of the chunk in dataset
        start_idx = i * chunk_size
        end_idx = (i + 1) * chunk_size

        # get the chunk of data
        redshift_ds_chunk = redshift_ds[start_idx:end_idx]
        BHpos_ds_chunk = BHpos_ds[start_idx:end_idx]
        BHMass_ds_chunk = BHMass_ds[start_idx:end_idx]

        # create z filter
        z_mask = np.abs(redshift_ds_chunk - zmerge) < 1e-2
        # create position filter
        target_BHpos_broadcast = np.broadcast_to(
            target_BHpos, (len(BHpos_ds_chunk[z_mask]), target_BHpos.shape[0]))
        pos_diff = BHpos_ds_chunk[z_mask] - target_BHpos_broadcast
        r = np.linalg.norm(pos_diff, axis=1)

        candidate_idx = np.where(
            (BHMass_ds_chunk[z_mask] == target_acBHMass) & (r < 5))[0]

        if len(candidate_idx) > 0:
            id2 = BHID_ds[start_idx:end_idx][z_mask][candidate_idx[0]]
            idx2 = start_idx + np.where(z_mask)[0][candidate_idx[0]]
            cur_chunk_i = i
            break

    if id2 is None:
        print(
            "Within dz=0.01, didn't find matched ID2 with acBHmass", flush=True)
    else:
        print(f"Found ID2 = {id2}, r = {r[candidate_idx[0]]:.2f}", flush=True)

    return id2, idx2



def validate_ID2_by_z(target_BHID, zmerge):
    """
    There are old swallowed ghost BHs still stay in BHdetailfile
    We need to check that ID2 is indeed the BH swallowed around zmerge
    """
    matching_indices = np.where(BHID_ds == target_BHID)
    target_z = redshift_ds[matching_indices]
    target_swallowed = swallowed_ds[matching_indices]

    z_mask2 = (target_z <= zmerge) & (target_z >= zend)
    criteria = (np.all(target_swallowed[target_z > zmerge] == 0)) & (
        np.all(target_swallowed[z_mask2] == 1))
    return criteria


def find_z_when_swallowed_flips(target_BHID):
    """
    Find the redshift where swallowed first flip 1 (for ID2)
    """
    matching_indices = (BHID_ds == target_BHID).nonzero()[0]

    target_z = redshift_ds[matching_indices]
    target_swallowed = swallowed_ds[matching_indices]

    sorted_indices = np.argsort(-target_z)
    target_z = target_z[sorted_indices]
    target_swallowed = target_swallowed[sorted_indices]

    first_positive_swallowed_index = np.argwhere(target_swallowed > 0)[0]
    idx2 = matching_indices[sorted_indices[first_positive_swallowed_index]]
    return target_z[first_positive_swallowed_index], idx2



def set_path():
    global snap, cluster
    if len(sys.argv) != 3:
        sys.exit("Please provide two arguments: cluster and snap number.")
    cluster = str(sys.argv[1]).lower()
    snap = str(sys.argv[2])

    if cluster == "vera":
        root = "/hildafs/datasets/Asterix"
        bh_reduce_file = f"{root}/BH_details_bigfile{'2' if snap > '347' else ''}/BH-Details-R{snap}"
        output_data_file = f"{root}/BH-mergers/bh-merger-extended-R{snap}"
    elif cluster == "frontera":
        bh_reduce_file = f"/home1/08942/nianyic/scratch3/Astrid/bhdetails/BH-Details-R{snap}"
        output_data_file = f"/home1/08942/nianyic/work/Astrid/mergers/extended/bh-merger-extended-R{snap}"
    elif cluster == "bridges":
        root = "/jet/home/nianyic/scratch1/Astrid/bhdetails-chopped"
        bh_reduce_file = f"{root}/BH-Details-R{snap}"
        output_data_file = f"{root}/mergers/bh-merger-extended-R{snap}"
    else:
        sys.exit("unknown cluster name. Please provide either vera or frontera.")
    return bh_reduce_file, output_data_file


def clean_memory():
    global swallowed_ds, swallowID_ds
    swallowed_ds = None
    swallowID_ds = None


def get_all_details_file():
    if cluster == "vera":
        root = "/hildafs/datasets/Asterix"
        flist = []
        dirlist = [f"{root}/BH_details_bigfile", f"{root}/BH_details_bigfile2"]
        for d in dirlist:
            flist = flist + sorted(glob.glob(f"{d}/BH-Details-R*"))
        flist.remove(f"{root}/BH_details_bigfile/BH-Details-R087")
        fidx_to_path = {i: f for i, f in enumerate(flist)}
        snap_to_fidx = {int(f.split("-R")[-1]): i for i, f in enumerate(flist)}
    elif cluster == "frontera":
        root = "/home1/08942/nianyic/scratch3/Astrid/bhdetails"                 
        flist = flist + sorted(glob.glob(f"{d}/BH-Details-R*"))            
        fidx_to_path = {i: f for i, f in enumerate(flist)}                     
        snap_to_fidx = {int(f.split("-R")[-1]): i for i, f in enumerate(flist)}
    elif cluster == "bridges":
        root = "/jet/home/nianyic/scratch1/Astrid/bhdetails-chopped"
        flist = sorted(glob.glob(f"{root}/BH-Details-R*"))
        fidx_to_path = {i: f for i, f in enumerate(flist)}
        snap_to_fidx = {f.split("-R")[-1]: i for i, f in enumerate(flist)}
    return fidx_to_path, snap_to_fidx


def main():
    bh_reduce_file, output_data_file = set_path()
    load_bh_data(bh_reduce_file)
    events = get_index_data()
    clean_memory()
    merger_data = get_bh_info(events)

    if merger_data is not None:
        np.save(output_data_file, merger_data)
        print(f"Saved merger data to {output_data_file}.npy", flush=True)


if __name__ == "__main__":
    main()
