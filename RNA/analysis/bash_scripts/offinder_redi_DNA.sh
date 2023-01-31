#!/bin/bash


source ../config.sh

while getopts "s:f:o:" flag
do
    case "${flag}" in
        s) sample=${OPTARG};;
        f) file=${OPTARG};;
        o) output=${OPTARG};;
    esac
done
python $reditools/src/cineca/reditools.py \
-B $file \
-f /g/steinmetz/project/dcm_WGS/wgs_variant_calling/bam/$sample.sort.dedup.recal.bam \
-r $reference_seq \
-o $output
