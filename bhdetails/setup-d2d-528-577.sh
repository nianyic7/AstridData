#!/bin/bash

for i in 528 534 539 540 541 547 552 557 561 566 572 577; do
    cp sample-d2d.sh job_d2d-$i.sh
    sed -i "s/placeholder/$i/g" job_d2d-$i.sh
done
