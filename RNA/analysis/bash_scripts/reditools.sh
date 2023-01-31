#!/bin/bash


source ../config.sh

while getopts s: flag
do
    case "${flag}" in
        s) sample=${OPTARG};;
    esac
done

python $reditools/src/cineca/reditools.py \
-f $basedir/bams/$sample.Aligned.sortedByCoord.out.dedup.split.recal.bam -g chr19:53843040-53843440 \
-r $basedir/genome/mm10_AAV.fa \
-o $basedir/bystander_edits/$sample.all_edits.txt


