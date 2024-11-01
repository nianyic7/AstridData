#!/bin/bash

for i in 632 638 641 642 647 648 653; do
    sbatch job_d2d-$i.sh
done
