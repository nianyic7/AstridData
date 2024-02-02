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

hh = 0.6774

acBHMass_ds = None
swallowed_ds = None
swallowID_ds = None
BHID_ds = None
redshift_ds = None
BHpos_ds = None
BHMdot_ds = None
BHMass_ds = None
zstart = None
zend = None


cur_chunk_i = 0




def load_bh_data(filename):
    global acBHMass_ds, swallowed_ds, swallowID_ds, BHID_ds, redshift_ds, BHpos_ds, BHMass_ds, zstart, zend

    print(f"Reading {filename}")
    bhd = BigFile(filename)
    # load partial data to check if there is merger event
    acBHMass_ds = bhd.open('acBHMass')[:]
    swallowed_ds = bhd.open('Swallowed')[:]
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
    BHpos_ds = bhd.open('BHpos')[:][sort_indices]
    BHMass_ds = bhd.open('BHMass')[:][sort_indices]
    BHMdot_ds = bhd.open('Mdot')[:][sort_indices]

    # get redshift range
    zstart, zend = np.max(redshift_ds), np.min(redshift_ds[redshift_ds > 0])
    print("load_bh_data, Total length:", len(BHID_ds))
    print("Redshift range:", zstart, zend, flush=True)



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
    ID2_to_index2, ID2_to_ID1 = get_swallowed_events() # indices of BH2 at merger, note: this may be a subset of ID1_to_index1

    idx1_ac = np.array(list(index1_to_ID1.keys()))
    ID1_sw = np.array(list(ID2_to_ID1.values()))
    ID2_sw = np.array(list(ID2_to_index2.keys()))

    num_mergers = 0
    data = []
    case1 = 0
    case2 = 0
    case3 = 0
    case4 = 0
    case5 = 0
    for idx1 in idx1_ac:
        id_1 = index1_to_ID1[idx1]
        assert id_1 == BHID_ds[idx1]

        zmerge = redshift_ds[idx1]
        # all events where ID1 is swallower
        swallowed_idx = np.where(ID1_sw == id_1)[0]

        id_2 = None
        if len(swallowed_idx) == 0: 
            # case1: BH1 swallowed the mass, but is not recoreded on BH2 history
            print(f"ID1={id_1}, Mass={BHMass_ds[idx1]}, acBHMass={acBHMass_ds[idx1]}, "
                  f"now trying to find ID2 by position and mass", flush=True)
            id_2, idx2 = find_ID2_by_pos_mass(zmerge, acBHMass_ds[idx1], BHpos_ds[idx1])
            case1 += 1

        elif len(swallowed_idx) == 1:
            # case2: the easy case, only one BH2 is swallowed and correctly recorded
            id_2 = ID2_sw[swallowed_idx[0]]
            idx2 = ID2_to_index2[id_2]
            case2 += 1

        else:
            # case3: multiple BH2 are swallowed, find the one logged and closest to z_acBH
            for potential_id in ID2_sw[swallowed_idx]:
                if validate_ID2_by_z(potential_id, zmerge):
                    id_2 = potential_id
                    idx2 = ID2_to_index2[id_2]
                    case3 += 1
                    break
        
        # if still no id_2, try to find it by zflip
        if ((len(swallowed_idx) > 1) & (id_2 is None)):
            zflip_idx2 = np.array([find_z_when_swallowed_flips(potential_id)
                            for potential_id in ID2_sw[swallowed_idx]])
            z_flip = zflip_idx2[:, 0]
            idx2 = zflip_idx2[:, 1]

            closest_z_flip_idx = np.argmin(np.abs(z_flip - zmerge))
            id_2 = ID2[swallowed_idx[closest_z_flip_idx]]
            idx2 = idx2[closest_z_flip_idx]
            print(
                f"With second zflip search, find ID2 = {id_2}, zflips = {z_flip}", flush=True)
            case4 += 1
        if id_2 is None:
            case5 += 1
            print(
                f"Still did not find matched a BH2, !!! acBHMass = {acBHMass_ds[i]}", flush=True)
            continue
        
        # convention is m1 > m2
        m1 = BHMass_ds[idx1] - acBHMass_ds[idx1]
        m2 = acBHMass_ds[idx1]
        if m1 < m2:
            m1, m2 = m2, m1
            id_1, id_2 = id_2, id_1
            idx1, idx2 = idx2, idx1
        data.append([id_1, id_2, idx1, idx2])
    print("Total number of mergers found:", len(data), flush=True)
    print(f"case1: {case1}, case2: {case2}, case3: {case3}, case4: {case4}, case5: {case5}", flush=True)


    #     merger_data[num_mergers] = (
    #         zmerge, id_1, id_2, m1 * 1e10 / hh, m2 * 1e10 / hh)
    #     num_mergers += 1

    # merger_data = merger_data[:num_mergers]
    # print(f"Check number of mergers: {len(merger_data)}")

    # merger_data = np.unique(merger_data)
    # print("Unique event number:", len(merger_data))

    return data


def get_bh_info(events):
    """_summary_

    Args:
        id1 (_type_): _description_
        id2 (_type_): _description_
        idx1 (_type_): _description_
        idx2 (_type_): _description_

    Returns:
        _type_: _description_
    """
        # initialize empty np array for merger_data

    merger_dtype = np.dtype(
        [('z', 'd'), ('ID1', 'q'), ('ID2', 'q'), ('m1', 'd'), ('m2', 'd'), ('mdot1', 'd'), \
            ('mdot2', 'd'), ('pos1', '3d'), ('pos2', '3d'), ('v1', '3d'), ('v2', '3d')])
    merger_data = np.empty(len(events), dtype=merger_dtype)

    for i, (id1, id2, idx1, idx2) in enumerate(events):
        m1 = BHMass_ds[idx1] - acBHMass_ds[idx1]
        m2 = acBHMass_ds[idx1]
        id1_1 = BHID_ds[idx1 - 1]
        id2_1 = BHID_ds[idx2 - 1]
        assert id1_1 == id1, "swallow happened at the very beginning of the chunk, cannot get Pos/Mdot info"
        assert id2_1 == id2, "swallow happened at the very beginning of the chunk, cannot get Pos/Mdot info"

        # validate mass info
        m1_alt = BHMass_ds[idx1 - 1]
        m2_alt = BHMass_ds[idx2 - 1]
        assert np.isclose(m1, m1_alt, rtol=1e-7), f"Mass info from before/after merger differs, m1 = {m1}, m1_alt = {m1_alt}"
        assert np.isclose(m2, m2_alt, rtol=1e-7), f"Mass info from before/after merger differs, m2 = {m2}, m2_alt = {m2_alt}"

        # get mdot info
        mdot1 = BHMdot_ds[idx1 - 1]
        mdot2 = BHMdot_ds[idx2 - 1]
        pos1 = BHpos_ds[idx1 - 1]
        pos2 = BHpos_ds[idx2 - 1]
        vel1 = BHvel_ds[idx1 - 1]
        vel2 = BHvel_ds[idx2 - 1]
        merger_data[i] = (zmerge, id_1, id_2, m1 * 1e10 / hh, m2 * 1e10 / hh, \
            mdot1, mdot2, pos1, pos2, vel1, vel2)

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
    if len(sys.argv) != 3:
        sys.exit("Please provide two arguments: cluster and snap number.")
    cluster = str(sys.argv[1]).lower()
    snap = str(sys.argv[2])

    if cluster == "vera":
        root = "/hildafs/datasets/Asterix"
        bh_reduce_file = f"{root}/BH_details_bigfile{'2' if snap > '347' else ''}/BH-Details-R{snap}"
        output_data_file = f"{root}/BH_mergers/bh-merger-extended-R{snap}"
    elif cluster == "frontera":
        bh_reduce_file = f"/scratch3/06431/yueyingn/BH-detail-reduce/BH-Details-R{snap}"
        output_data_file = f"mergers_below2/bh-merger-extended-R{snap}"
    else:
        sys.exit("unknown cluster name. Please provide either vera or frontera.")
    return bh_reduce_file, output_data_file




def main():
    bh_reduce_file, output_data_file = set_path()
    load_bh_data(bh_reduce_file)
    events = get_index_data()
    merger_data = get_bh_info(events)



    # if merger_data is not None:
    #     np.save(output_data_file, merger_data)


if __name__ == "__main__":
    main()
