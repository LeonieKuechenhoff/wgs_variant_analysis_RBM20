#!/bin/bash

#conda activate mamba

#script to run reditools on potential off target positions with DNA data

source ../config.sh

#DNA
for sample in L279 L450 L282 H279 H450 H282 T279 T450 T282
do
    sbatch --job-name=pspry -t 7-00:00:00 --cpus-per-task=12 --mem=150g --output=off_finder_redi_out/redi_DNA_$sample.out \
    ./offinder_redi_DNA.sh -s $sample -f $basedir/offinder/off_target_spry.bed \
    -o $basedir/offinder/$sample.DNA.redi.txt
done
