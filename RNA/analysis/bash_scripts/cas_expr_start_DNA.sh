#!/bin/bash

module load deepTools/3.5.0-foss-2021a
source ../config.sh


for sample in H279 L279 T279 H282 L282 T282 H450 L450 T450
do
    sbatch --job-name=cov -t 3-00:00:00 --cpus-per-task=4 --mem=10g --output=$basedir/jobs/cas_expr/$sample.dna.out \
    ./cas_exprDNA.sh -s $sample
done

