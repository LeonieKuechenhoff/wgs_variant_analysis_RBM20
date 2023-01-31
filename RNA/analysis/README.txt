Files

- config.py 
Paths that are used several times in scripts. Are imported at start of most scripts
- config.sh
Paths that are used several times in bash scripts.

- tissue_spec_list.ipynb
Write one list per tissue of tissue specific variants (with anntation in which other samples these variant occur)

- filter_mgp_variants.ipynb
Filter mgp db entries based on their annotations. Variants of interest are only protein coding variants.
- overlap_callers.ipynb
Compare outputs of 6 different varaint callers. Determine count of called variants, overlap with protein coding known variants (extracted with script above)
an overlap with each other

- variant_caller_venn_HC_PL_ST.ipynb
Script to plot venn diagrams for variants called by HC, PL and ST 
-on_target_AF.pynb
Plot allele frequency of on-target edits
- AF.ipynb 
Plot counts of tissue specific and common variants & check relationship to read coverage
-SNP_type_pm_strand.ipynb
Calculate fraction of SNP per sample and categorize into type of SNP

-annotation.ipynb
 Annotate txt files and add info if variant is on pos. or negative strand. Save annotated df in annotation_dir

 - bystander_edits.ipynb
 Search for bystander edits (edits in gRNA binding region)

 - bystander_edits_DNA.ipynb
 Search for bystander edits (edits in gRNA binding region) but with DNA data

- variant_sourrounding.ipynbScript to check list of variats in more detail.
  The following things will be checked:
Are there any gRNA Sequence similarities in proximity?

- cas_offinder.ipynb
Import potential off-target sites (found by cas-offinder) & export bam file with regions of interest

- offfinder_redi.ipynb 
Import reads exported with reditools in all potential off-target sites (found by cas-offinder, <= 5 mismatches to gRNA) and plot to see a trend

- compare_cov.ipynb
Script to compare coverge across samples




Bash scripts
Syntax:
a) Script to start script b on slurm
b) Actual script


a)coverage_start.sh
b)coverage.sh
script to extract read coverage in 100 bp bins 

a) offinder_start.sh
b) offinder.sh
script to extract potentail off-target sites by getting sites with few mismatches across the genome

a)offinder_redi_start_DNA.sh
b)offinder_redi_DNA.sh
script to run reditools on potential off target positions with RNA data

a) offinder_redi_start.sh
b) offinder_redi.sh
script to run reditools on potential off target positions with DNA data


a) reditools_start.sh
b) reditools.sh
script to run reditools on area around target-site with RNA data

readcount.sh  
Script to give read count of primary aligned mapped reads per sample
