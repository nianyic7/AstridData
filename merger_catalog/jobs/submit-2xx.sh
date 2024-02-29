#!/bin/bash
for f in job_merger_2*.sh
do 
    echo "Processing $f file..."
    sbatch $f
    sleep 2
done
