"""
extract merger event from bhdetail-reduce files
Usage: python extract-merger-v1.py $snap

This script is used to extract merger events from BHdetail in BigFile formatting.
It will output a npy file containing merger data,
for each merger event, it contains [redshift, ID1, ID2, m1, m2].
"""

import numpy as np
from bigfile import BigFile
import sys

hh = 0.6774

acBHMass_ds = None
swallowed_ds = None
swallowedID_ds = None
BHID_ds = None
redshift_ds = None
BHpos_ds = None
BHMass_ds = None
zstart = None
zend = None


cur_chunk_i = 0


def get_merger_data(filename):
    global acBHMass_ds, swallowed_ds, swallowedID_ds, BHID_ds, redshift_ds, BHpos_ds, BHMass_ds, zstart, zend

    bhd = BigFile(filename)

    # load partial data to check if there is merger event
    acBHMass_ds = bhd.open('acBHMass')[:]
    swallowed_ds = bhd.open('Swallowed')[:]
    # check if no merger event
    mask0 = (acBHMass_ds > 1e-8) & (swallowed_ds == 0)
    if np.sum(mask0) < 1:
        print("No merger event here")
        return None

    # load full data
    swallowID_ds = bhd.open('SwallowID')[:]
    BHID_ds = bhd.open('BHID')[:]
    redshift_ds = bhd.open('z')[:]
    BHpos_ds = bhd.open('BHpos')[:]
    BHMass_ds = bhd.open('BHMass')[:]

    # get redshift range
    zstart, zend = np.max(redshift_ds), np.min(redshift_ds[redshift_ds > 0])
    print("Redshift range:", zstart, zend)

    # get ids
    mask1 = swallowID_ds > 0
    ID1 = swallowID_ds[mask1]  # swallower
    ID2 = BHID_ds[mask1]  # the one got swallowed

    # one BH can not be swallowed twice, so reduce ID2
    ID2, i = np.unique(ID2, return_index=True)
    ID1 = ID1[i]

    print(f"Number of mergers from acBHMass: {len(acBHMass_ds[mask0])}")
    print("Find unique SwallowID-BHID pairs", len(ID1))
    print("Number of unique SwallowID", len(np.unique(ID1)), flush=True)

    # start extracting merger data
    valid_indices = np.where(mask0)[0]

    # initialize empty np array for merger_data
    merger_dtype = np.dtype(
        [('z', 'd'), ('ID1', 'q'), ('ID2', 'q'), ('m1', 'd'), ('m2', 'd')])
    merger_data = np.empty(len(valid_indices), dtype=merger_dtype)

    num_mergers = 0

    for i in valid_indices:
        zmerge = redshift_ds[i]
        id_1 = BHID_ds[i]  # swallower
        swallowed_idx = np.where(ID1 == id_1)[0]

        id_2 = None

        if len(swallowed_idx) == 1:
            id_2 = ID2[swallowed_idx[0]]
        elif len(swallowed_idx) < 1:
            print(f"No matched swallowID, ID1={id_1}, Mass={bhd['BHMass'][i]}, acBHMass={acBHMass_ds[i]}, "
                  f"now trying to find ID2 by position and mass", flush=True)
            id_2 = find_ID2_by_pos_mass(zmerge, acBHMass_ds[i], BHpos_ds[i])
        else:
            # multiple ID2, find the one that matches z
            for potential_id in ID2[swallowed_idx]:
                if validate_ID2_by_z(potential_id, zmerge):
                    id_2 = potential_id
                    # break
                    # TODO: og no break here, why?

        # if still no id_2, try to find it by zflip
        if ((len(swallowed_idx) > 1) & (id_2 is None)):
            z_flip = np.array([find_z_when_swallowed_flips(potential_id)
                               for potential_id in ID2[swallowed_idx]])
            closest_z_flip_idx = np.argmin(np.abs(z_flip - zmerge))
            id_2 = ID2[swallowed_idx[closest_z_flip_idx]]
            print(
                f"With second zflip search, find ID2 = {id_2}, zflips = {z_flip}", flush=True)

        if id_2 is None:
            print(
                f"Still did not find matched a BH2, !!! acBHMass = {acBHMass_ds[i]}", flush=True)
            continue

        m1 = BHMass_ds[i] - acBHMass_ds[i]
        m2 = acBHMass_ds[i]
        # convention is m1 > m2
        if m1 < m2:
            m1, m2 = m2, m1
            id_1, id_2 = id_2, id_1

        merger_data[num_mergers] = (
            zmerge, id_1, id_2, m1 * 1e10 / hh, m2 * 1e10 / hh)
        num_mergers += 1

    merger_data = merger_data[:num_mergers]
    print(f"Check number of mergers: {len(merger_data)}")

    merger_data = np.unique(merger_data)
    print("Unique event number:", len(merger_data))

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
            cur_chunk_i = i
            break

    if id2 is None:
        print(
            "Within dz=0.01, didn't find matched ID2 with acBHmass", flush=True)
    else:
        print(f"Found ID2 = {id2}, r = {r[candidate_idx[0]]:.2f}", flush=True)

    return id2


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
    matching_indices = np.where(BHID_ds == target_BHID)

    target_z = redshift_ds[matching_indices]
    target_swallowed = swallowed_ds[matching_indices]

    sorted_indices = np.argsort(-target_z)
    target_z = target_z[sorted_indices]
    target_swallowed = target_swallowed[sorted_indices]

    first_positive_swallowed_index = np.argwhere(target_swallowed > 0)[0]
    return target_z[first_positive_swallowed_index]


def main():
    if len(sys.argv) != 2:
        print("Need to specify snapshot number")
        return

    snap = str(sys.argv[1])
    # TODO: path should be an argument?
    bh_reduce_file = f"/hildafs/datasets/Asterix/BH_details_bigfile{'2' if snap > '347' else ''}/BH-Details-R{snap}"
    output_data_file = f"mergers_below2/bh-merger-R{snap}"

    print(f"Reading {bh_reduce_file}")

    merger_data = get_merger_data(bh_reduce_file)

    if merger_data is not None:
        np.save(output_data_file, merger_data)


if __name__ == "__main__":
    main()
