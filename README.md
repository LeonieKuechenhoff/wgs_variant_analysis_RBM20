# wgs_variant_analysis_RBM20
Whole genome sequencing analysis to investigate on- and off-target effects of CRISPR base-editor treatment.

## Folder structure
### 1. /variant_calling
Pre-processing of the raw fastq files into called variants (.vcf).  
#### Steps:
##### 1. Install and activate conda environment
mamba env create -p ./envs/smake -f ./envs/environment.yaml  
conda activate ./envs/smake
##### 2. Prepare input and configuration files
inputs/input.yaml  
inputs/samples.txt  
inputs/fasta/mm10_AAV.fa.gz ### reference genome  
inputs/fasta/mm10_AAV.bed ### bed file containing chromosomal lengths required for scalpel  
inputs/vcf/*.vcf.gz ### known sites for BQSR,VQSR and lofreq  
##### 3. Run snakemake
snakemake --profile ./profile  

### 2. /variantfile_preparation
Steps following variant calling. Includes snakemake workflow to normalize, filter, and merge vcf files from differnt tissues and variant callers into one txt file.   
a) snakefile  
    - Snakemake workflow  
b) /scripts  
    i) vcf_to_AFtable.pl  
        Transforms vcf files to txt table for easier data handling.  
    ii) db_entry_overlap.py  
        Script to extract overlapping sites with sites annotated in databases. Also corrects genotype of duplicated alleles.  
Returns trimmed database entry files with only entries overlapping with sample variants.  
    iii) create_filter_tables.py  
        Script to perform quality filtering and save filtered variants  
    iv) table_annovar.pl  
        Script to add genome annotation to variants.  

### 3. /variant_analysis
Analysis of serveral aspects of output files from step two.
config.py - directory settings for import statements. Most import files are created in step 2.  
a) variant_caller_venn.ipynb  
    Plot venn diagrams of variants called by different variant callers.  
b) tissue_spec_list.ipynb  
    Write one list per tissue of tissue specific variants. Output needed for AF.ipynb & variant_surrounding.ipynb  
c) AF.ipynb  
    Plot total number of tissue specific and common variants  
    Plot allele frequency of A to G mutations  
    Plot allele frequency of on-target edits  
d) annotation.ipynb  
    Plot how variants are annotated (intron vs. exon etc.), caluclate fraction per category, sample and tissue.  
e) SNP_type.ipynb  
    Calculate fraction of SNP per sample and tissue and categorize into type of SNP, plot results.  
f) variant_surrounding.ipynb  
    For each tissue specific variant, the following things will be checked:  
    Is it an A>G mutation?  
    Are there any gRNA Sequence similarities around the variant?  
g) semi_global_alignment.pl  
    Helper scrip that is called in variant_surrounding.ipynb. Searches for Seqeuence similarities of two sequences and  
    returns number of mismatches of best match.
