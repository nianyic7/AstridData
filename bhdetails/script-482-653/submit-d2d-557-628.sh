#!/bin/bash

for i in 557 561 566 572 577 582 583 586 592 594 599 600 604 606 612 616 621 622 627 628; do
    sbatch job_d2d-$i.sh
done
