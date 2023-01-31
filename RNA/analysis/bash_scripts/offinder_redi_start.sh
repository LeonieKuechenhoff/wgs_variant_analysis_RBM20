#!/bin/bash

#conda activate mamba

#script to run reditools on potential off target positions with DNA data

source ../config.sh

#RNA
for sample in H029_pbs_R L029_pbs_R H028_pbs_R L028_pbs_R H032_pbs_R L032_pbs_R H033_nrch_R L033_nrch_R H030_nrch_R L030_nrch_R L036_nrch_R H036_nrch_R
do
    sbatch --job-name=rnrch -t 7-00:00:00 --cpus-per-task=8 --mem=50g --output=off_finder_redi_out/redi_rnrch_$sample.out \
    ./offinder_redi.sh -s $sample -f $basedir/offinder/off_target_nrch_r.bed \
    -o $basedir/offinder/$sample.rnrch.redi.txt
done


for sample in L012_nrch L013_nrch L014_nrch H012_nrch H013_nrch H014_nrch L011_pbs L321_pbs L333_pbs H011_pbs H321_pbs H333_pbs
do
    sbatch --job-name=pnrch -t 7-00:00:00 --cpus-per-task=8 --mem=50g --output=off_finder_redi_out/redi_pnrch_$sample.out \
    ./offinder_redi.sh -s $sample -f $basedir/offinder/off_target_nrch_r.bed \
    -o $basedir/offinder/$sample.pnrch.redi.txt
done


for sample in L279_spry L450_spry L283_spry H279_spry H283_spry H450_spry L011_pbs L321_pbs L333_pbs H011_pbs H321_pbs H333_pbs
do
    sbatch --job-name=pspry -t 7-00:00:00 --cpus-per-task=8 --mem=50g --output=off_finder_redi_out/redi_pspry_$sample.out \
    ./offinder_redi.sh -s $sample -f $basedir/offinder/off_target_spry.bed \
    -o $basedir/offinder/$sample.pspry.redi.txt
done

