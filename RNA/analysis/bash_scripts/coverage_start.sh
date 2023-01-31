#!/bin/bash
#module load deepTools/3.5.0-foss-2021a

for sample in H029_pbs_R L029_pbs_R H028_pbs_R L028_pbs_R H032_pbs_R L032_pbs_R H033_nrch_R L033_nrch_R H030_nrch_R L030_nrch_R L036_nrch_R H036_nrch_R L012_nrch L013_nrch L014_nrch H012_nrch H013_nrch H014_nrch L279_spry L450_spry L283_spry H279_spry H283_spry H450_spry L011_pbs L321_pbs L333_pbs H011_pbs H321_pbs H333_pbs
do
    sbatch --job-name=cov -t 3-00:00:00 --cpus-per-task=8 --mem=100g --output=cov_out/$sample.out \
    ./coverage.sh -s $sample
done

