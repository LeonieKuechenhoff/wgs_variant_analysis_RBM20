#!/bin/bash

# script to extract read coverage in 100 bp bins 
source ../config.sh

while getopts "s:" flag
do
    case "${flag}" in
        s) sample=${OPTARG};;
    esac
done
bamCoverage --bam /g/steinmetz/project/dcm_WGS/analysis/bams/$sample.sort.dedup.recal.bam \
--outFileName $basedir/cas_expr/$sample.DNA.bed --outFileFormat bedgraph --normalizeUsing CPM --region AAV_Nterm_SpRY:5126:6841 --binSize 1