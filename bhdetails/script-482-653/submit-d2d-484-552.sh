#!/bin/bash

for i in 484 486 491 493 496 497 500 506 508 509 514 518 523 528 534 539 540 541 547 552; do
    sbatch job_d2d-$i.sh
done
