#!/bin/bash

for i in 528 534 539 540 541 547 552 557 561 566 572 577; do
    cp d2d-sample.sh d2d-$i.sh
    sed -i "s/placeholder/$i/g" d2d-$i.sh
done
