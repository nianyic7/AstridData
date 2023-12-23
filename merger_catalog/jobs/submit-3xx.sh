#!/bin/bash
for f in job_merger_3*.sh
do 
    echo "Processing $f file..."
    sbatch $f
    sleep 2
done
