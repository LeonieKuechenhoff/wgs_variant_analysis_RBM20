#!/bin/bash
source ../config.sh

for sample in H011_pbs H012_nrch H013_nrch H014_nrch H279_spry H283_spry H321_pbs H333_pbs H450_spry H029_pbs_R L029_pbs_R H033_nrch_R L033_nrch_R L011_pbs L012_nrch L013_nrch L014_nrch L279_spry L321_pbs L333_pbs L450_spry L283_spry H028_pbs_R H030_nrch_R H032_pbs_R H036_nrch_R L028_pbs_R L030_nrch_R L032_pbs_R L036_nrch_R
do
    sbatch --job-name=strelka -t 7-00:00:00 --cpus-per-task=8 --mem=20g --output=$basedir/variant_caller_out/vartest/redi$sample.out \
    ./reditools.sh -s $sample
done




