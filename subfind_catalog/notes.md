## ASTRID Subgroup Note
**Chunk data**
- group in chunks:
    chunki.j exclusive of groupidx j (max pig group idx in this chunk is j-1!)

- Are some PIG groups not in Subfind groups?
    Yes but not many (see `subfind_test.py`), ~ 1e3 for large groups.
    
- So what do we do about them?
   only follow the original FOF halo, discard the subfind groups and use subgroups only. Hopefully (checked in the code and did not find any so far) no subgroups are splitted into two FOF groups.


## Pipeline
`save_subfind_order.py`: save a chunkwise-unique subhalo ID for each particle in the new PIG file, following the subID-sorted order. Save a `NewIndex` column that maps the original Groups' particle ordering into the subID-sorted order. Can use `np.argsort(NewIndex)` to revert the mapping.

`group_reorder.py`: For all columns in the original PIG file, write the re-ordered version for each group following the "NewOrder" column; 

**NOTE:** the second half of the groups contains a small amount of particles, but takes up a huge chunk of computation time; to save time, we do not reorder anything when the total particle in group is below 100 (i.e. we assume that those halos do not contain subhalos, which should be reasonable)!


`subgroup_summary.py`: compute/write the  `SubGroups/LengthByType` and `SubGroups/OffsetByType` column; compute/write the  `FOFGroups/GroupNSub` and `FOFGroups/FirstSub` column; migrate all subgroup properties from chunks in the new ordering.

## Caution
- Can work on a subset of chunks, but remember NOT to create existing blocks again, or the existing data may get erased. Just to `dest[blockname].write(start,data)`.
- Later chunks (second half) are much slower to process because each contains too many groups. Allocate more processors to the second half to balance the work.
 
 

## Work Log
Oct 11

- PIG_214 (z=3) all processed and well-tested

Oct 18

- Fixed a bug in `group_reorder.py`: assumed all groups has DM particle which is not true. Affect chunk>810, need to reprocess those, should be quick!
 
 Oct 19
- No more bugs in the chunk->BF code.
- However for PIG_348 (z=2): post-processing code works fine, but the bug is within Subfind itself (e.g. chunk0.1 after subgroup 232). Need to resolve this issue before processing again.