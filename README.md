# wgs_variant_analysis_RBM20
Whole genome and transcriptome sequencing analysis to investigate on- and off-target effects of CRISPR base-editor treatment.

Publication: 
Grosch, M., Schraft, L., Chan, A. et al.  
Striated muscle-specific base editing enables correction of mutations causing dilated cardiomyopathy.   
Nat Commun 14, 3714 (2023).   
https://doi.org/10.1038/s41467-023-39352-1


## DNA

### 1. /variant_calling
Pre-processing of the raw fastq files into called variants (.vcf).  
#### Steps:
##### 1. Install and activate conda environment
mamba env create -p ./envs/smake -f ./envs/environment.yaml  
conda activate ./envs/smake  
##### 2. Prepare input and configuration files
./variant_calling/RunThisFirst.sh
##### 3. Run snakemake
snakemake --profile ./profile  

### 2. /variantfile_preparation
Steps following variant calling. Includes snakemake workflow to normalize, filter, and merge vcf files from differnt tissues and variant callers into one txt file. Can be run from the same environment as the variant calling part. Therefore, call:   
conda activate ./envs/smake

#### a) snakefile  
Contains snakemake workflow. Call with:  
snakemake --profile profile
#### b) /scripts  
Helper scripts that are called in the snakemake workflow:  
 - vcf_to_AFtable.pl  
Transforms vcf files to txt table for easier data handling.    
 - db_entry_overlap.py  
Script to extract overlapping sites with sites annotated in databases. Also corrects genotype of duplicated alleles.  
Returns trimmed database entry files with only entries overlapping with sample variants.  
 - create_filter_tables.py  
Script to perform quality filtering and save filtered variants  
 - table_annovar.pl  
Script to add genome annotation to variants.  

### 3. /variant_analysis
Analysis of serveral aspects of output files from step two.
config.py - directory settings for import statements. Most import files are created in step 2.  
 - variant_caller_venn.ipynb 
    Plot venn diagrams of variants called by different variant callers.  
 - tissue_spec_list.ipynb  
    Write one list per tissue of tissue specific variants. Output needed for AF.ipynb & variant_surrounding.ipynb  
 - AF.ipynb  
    Plot total number of tissue specific and common variants  
    Plot allele frequency of A to G mutations  
    Plot allele frequency of on-target edits   
 - annotation.ipynb    
    Plot how variants are annotated (intron vs. exon etc.), calculate fraction per category, sample and tissue.  
 - SNP_type.ipynb  
    Calculate fraction of SNP per sample and tissue and categorize into type of SNP, plot results.  
 - variant_surrounding.ipynb    
    For each tissue specific variant, the following things will be checked:  
    Is it an A>G mutation?  
    Are there any gRNA Sequence similarities around the variant?  
 - semi_global_alignment.pl  
    Helper scrip that is called in variant_surrounding.ipynb. Searches for Seqeuence similarities of two sequences and  
    returns number of mismatches of best match.

## RNA
### 1. ./RNA/calling_prep
Pre-processing of the raw fastq files into variant files (.vcf). Also includes some steps following variant calling, such as filtering, normalization and merging them into a common txt file.
Can be run from the same environment as the DNA variant calling part. However certain rules of the snakemake pipeline require python 2.7, which is run from a separate environment. 
##### 1. Install python 2.7 conda environment. 
mamba env create -p ./env/smakep27 -f ./env/env_p27.yml  
The pipeline is run from the same evironment as the DNA pipeline. The p27 environment will be called withing the snakemake pipeline. Therefore, run:  
conda activate ../../DNA/variant_calling/env/smake
##### 2. Clone and follow installation instructions from the following git repositories:
 - opossum from https://github.com/BSGOxford/Opossum  
Add file location to ../../variantfile_preparation/config/config.json  
 - platypus from https://github.com/andyrimmer/Platypus  
Add file location to ../../variantfile_preparation/config/config.json  
 - strelka from https://github.com/Illumina/strelka  
Add file location to ../../variantfile_preparation/config/config.json  


##### 3. Run the snakemake pipeline
If all samples with the P635L mutation should be run, call:  
snakemake -s snakefile --profile profile --use-conda --config mutation=p635l  
If all samples with the R636Q mutation should be run, call:  
snakemake -s snakefile --profile profile --use-conda --config mutation=r636q  

### analysis  
This folder cotains all costum made scripts to analyze variants.  

##### 2. Clone and follow installation instructions from the following git repositories:
 - reditools2 from https://github.com/BioinfoUNIBA/REDItools2
Add file location to RNA/analysis/config.sh

The analysis of the data includes the following files:

#### Python scripts and jupyter notebooks

- config.py   
Paths that are used several times in scripts. Are imported at start of most scripts
- config.sh 
Paths that are used several times in bash scripts.  

- annotation.ipynb  
 Annotate txt files and add info if variant is on pos. or negative strand. Save annotated df in annotation_dir

- tissue_spec_list.ipynb 
Write one list per tissue of tissue specific variants (with anntation in which other samples these variant occur)

- AF.ipynb  
Plot counts of tissue specific and common variants & check relationship to read coverage  

- SNP_type_pm_strand.ipynb  
Calculate fraction of SNP per sample and categorize into type of SNP

 - bystander_edits.ipynb  
 Search for bystander edits (edits in gRNA binding region)

- variant_sourrounding.ipynb  
Script to check list of variants in more detail.  
Question that is adressed:  
Are there any gRNA Sequence similarities in proximity of each heart specific variant?

- semi_global_alignment.pl  
Helper scrip that is called in variant_surrounding.ipynb. Searches for Seqeuence similarities of two sequences and  
returns number of mismatches of best match.

- cass9_expr.ipynb
Script to compare nCas9 expression (requires cas_expr_start.sh to run first)

- signif_testing.ipynb
Script to test significance of AG variant proportion




#### Bash scripts - /analysis/bash/
Syntax:  
a) Script to start script b on slurm  
b) Actual script  


a)cas_expr_start.sh
b)cas_expr.sh
script to identify Cas9 expression

a) reditools_start.sh  
b) reditools.sh  
script to run reditools on area around target-site with RNA data  

readcount.sh
Script to give read count of primary aligned mapped reads per sample

