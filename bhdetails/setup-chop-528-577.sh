#!/bin/bash

for i in 528 534 539 540 541 547 552 557 561 566 572 577; do
    cp sample-chop.sh job_chop-$i.sh
    sed -i "s/placeholder/$i/g" job_chop-$i.sh
done
