#!/bin/bash
snap=483
dataroot="/home1/08942/nianyic/scratch3/Astrid"
pigfile="/home1/08942/nianyic/ASTRID2_PIG/PIG_$snap"
dest="$dataroot/PIG_483_subfind"
echo "copying header"
cp -r $pigfile/Header $dest
echo "copying fofgroups"
cp -r $pigfile/FOFGroups/* $dest/FOFGroups
echo "copying 4..."
cp -r $pigfile/4/GroupID $dest/4
echo "copying 5..."
cp -r $pigfile/5/GroupID $dest/5
echo "copying 0..."
cp -r $pigfile/0/GroupID $dest/0
echo "copying 1..."
cp -r $pigfile/1/GroupID $dest/1

echo "done!"
