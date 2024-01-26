def get_subLength(gidx,SidxFOF2Chunk,subLengthChunk):
    """
    Get the sorting index to sort by subgroup IDs
    for particles in each 
    """
    
    start,end = Offset[gidx], Offset[gidx] + Length[gidx]
    # assume that each subgroup must have a DM particle
    p = 1
    subID = dest_r['%d/SubgroupID'%p][start[p]:end[p]]
    subID_set   = set(subID)
    if len(subID_set) == 0 or min(list(subID_set)) == FUZZ:
        # no substructure in this group
        empty = np.zeros((0,6),dtype=np.int32)
        return np.array([],dtype=np.int32), empty, empty
    
    if max(list(subID_set)) == FUZZ:
        # there is inner fuzz
        subID_set.remove(FUZZ)
    subID_arr = np.sort(np.array(list(subID_set)))

    # N subgroups in this fof

        
    nsubs     = len(subID_set)
    subLength = np.zeros((nsubs,6),dtype='i4')
    subOffset = np.zeros((nsubs,6),dtype='i4')
    
    if len(subLengthChunk) == 0:
        assert nsubs == 0, 'chunk has no subgroup but group has! Gidx %d'%gidx
    
    for p in [0,1,4,5]:
        send = start[p]
        subID = dest_r['%d/SubgroupID'%p][start[p]:end[p]]
        for i,s in enumerate(subID_arr):
            mask,   = np.where(subID == s)
            slength = len(mask)
            
            # check alignment with the original subgroup
            assert slength == subLengthChunk[SidxFOF2Chunk[s]][p], \
            "subfind part number %d does not match with fof subgroup %d, gidx%d"\
            %(subLengthChunk[SidxFOF2Chunk[s]][p],slength,gidx)
            
            sstart = send
            send   = sstart + slength
            
            subLength[i][p] = slength
            subOffset[i][p] = sstart      
    return subID_arr, subLength, subOffset




def process_chunk(c):
    """
    For each chunk do the following:
        1. write all particle data of groups in new order
        2. obtain the lbt/obt for subgroups
        3. obtain the mapping from fof subgroup idx to chunk subgroups
            (just need to sort by mass)
            
    """
    # groups in this chunk
    print('Processing chunk:',c,flush=True)
    subdir  = subroot + '/chunk%d.%d/output/'%(c,maxgroup_list[c])
    tabfile = subdir + tab
    grpfile = subdir + grp

    if c == 0:
        gstart, gend = 0, maxgroup_list[c]
    else:
        gstart, gend = maxgroup_list[c-1], maxgroup_list[c]
        
    if "Subhalo" in h5py.File(tabfile,'r').keys():
        try:
            subsChunk  = h5py.File(tabfile,'r')['Header'].attrs['Nsubgroups_Total']
        except KeyError:
            subsChunk  = h5py.File(tabfile,'r')['Header'].attrs['Nsubhalos_Total']
            
        subLengthChunk = h5py.File(tabfile,'r')['Subhalo']['SubhaloLenType'][:]
        SubhaloTotMass = h5py.File(tabfile,'r')['Subhalo']['SubhaloMass'][:]
    else:
        subsChunk      = 0
        subLengthChunk = np.zeros((0,6),dtype=np.int32)
        SubhaloTotMass = np.zeros(0,dtype=np.float64)

        
        
    SidxFOF2Chunk = np.argsort(SubhaloTotMass)[::-1]
    
    NsubTot = 0
    subIDsetAll = set()
    subIDAll = []
    subLengthAll = []
    subOffsetAll = []
    for gidx in range(gstart,gend):
        # write features
        write_group(gidx)
        # get L/O
        subID_arr, subLength, subOffset = get_subLength(gidx,SidxFOF2Chunk,subLengthChunk)
        subID_set = set(subID_arr)
        if subIDsetAll & subID_set:
            raise Exception("two FOF halos contain particles from the same subgroup %d"%gidx)

        subIDsetAll.update(subID_set)
        subIDAll.append(subID_arr)
        subLengthAll.append(subLength)
        subOffsetAll.append(subOffset)
        
    subIDAll     = np.concatenate(subIDAll)
    subLengthAll = np.concatenate(subLengthAll,axis=0).astype(np.int32)
    subOffsetAll = np.concatenate(subOffsetAll,axis=0).astype(np.int32)
    ChunkSubID   = SidxFOF2Chunk[subIDAll]
    
    return subIDAll,subLengthAll,subOffsetAll,ChunkSubID


def init_all_blocks(source,dest,blocknames):
    pig = BigFile(source)
    blocklists = {p:[] for p in [0,1,4,5]}
    for p in [0,1,4,5]:
        features = get_features(pigdir=source,ptype=p)
        for f in features:
            blockname = '%d/%s'%(p,f)
            dtype,dsize,nfile = dtype_size_nfile(pig,blockname)
            block = dest.create(blockname, dtype, dsize, nfile)
            if rank == 0:
                print('Initizlized block:', blockname,flush=True)
            blocklists[p].append((block, blockname))
    return blocklists


                if len(data) > buffer_size:
                    if last == 0:
                    # buffer empty, directly write data
                        block.write(start[p],data)
                        buffer_start = end[p]
                    else:
                        # flush the buffer first
                        block.write(buffer_start, buffer[:last])
                        last = 0
                        block.write(start[p],data)
                        buffer_start = end[p]
                    
                elif len(data) <= buffer_size - last:
                    # append data to buffer
                    buffer[last:last+len(data)] = data
                    last += len(data)

                else:
                    # flush then add data to buffer
                    block.write(buffer_start,buffer[:last])
                    buffer[:len(data)] = data
                    buffer_start = start[p]
                    last = len(data)