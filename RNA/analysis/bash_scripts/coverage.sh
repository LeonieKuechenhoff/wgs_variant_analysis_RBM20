#!/bin/bash

# script to extract read coverage in 100 bp bins 
source ../config.sh

while getopts "s:" flag
do
    case "${flag}" in
        s) sample=${OPTARG};;
    esac
done
bamCoverage --bam $basedir/bams/$sample.Aligned.sortedByCoord.out.dedup.split.recal.bam \
--outFileName $basedir/coverage/$sample.bed --outFileFormat bedgraph