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
-f $basedir/bams/$sample.Aligned.sortedByCoord.out.dedup.split.recal.bam \
-r $basedir/genome/mm10_AAV.fa \
-o $output
