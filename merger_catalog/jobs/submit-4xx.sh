#!/bin/bash
for f in job_merger_4*.sh
do 
    echo "Processing $f file..."
    sbatch $f
    sleep 2
done
