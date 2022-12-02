#### Run this script to get external datasets required for reproducing the results for mice variant calling ####
### downloading mm10 reference genome ###
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz -O inputs/fasta/mm10.fa.gz
### combine mm10 reference with AAV plasmid sequence ###
zcat inputs/fasta/mm10.fa.gz inputs/fasta/AAV_ref.fa.gz > inputs/fasta/mm10_AAV.fa
bgzip inputs/fasta/mm10_AAV.fa
gatk CreateSequenceDictionary -R inputs/fasta/mm10_AAV.fa.gz

### get gtf annotation ###
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.refGene.gtf.gz -O inputs/fasta/mm10.refGene.gtf.gz

### get known SNPs ###
wget ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/00-All.vcf.gz -O inputs/vcf/mouse.dbsnp.vcf.gz
zcat inputs/vcf/mouse.dbsnp.vcf.gz | head -1000 | grep "^#" | cut -f 1-8 | sed 's/ID=\([0-9]\+\),/ID=chr\1,/' | sed 's/ID=\([XY]\),/ID=chr\1,/' | sed 's/ID=MT,/ID=chrM,/' > inputs/vcf/mouse.dbsnp.chr.vcf
zcat inputs/vcf/mouse.dbsnp.vcf.gz | grep -v "^#" | cut -f 1-8 | sed 's/^/chr/' | sed 's/chrMT/chrM/' >> inputs/vcf/mouse.dbsnp.chr.vcf
gatk UpdateVcfSequenceDictionary -SD inputs/vcf/mm10_AAV.dict -I inputs/vcf/mouse.dbsnp.chr.vcf -O inputs/vcf/mouse.dbsnp.chr.update.vcf
gatk SortVcf -SD inputs/vcf/mm10_AAV.dict -I inputs/vcf/mouse.dbsnp.chr.update.vcf -O inputs/vcf/mouse.dbsnp.chr.sort.vcf
bgzip inputs/vcf/mouse.dbsnp.chr.sort.vcf
tabix inputs/vcf/mouse.dbsnp.chr.sort.vcf.gz
rm inputs/vcf/mouse.dbsnp.chr.vcf inputs/vcf/mouse.dbsnp.chr.update.vcf

wget ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz -O inputs/vcf/mgp.v5.merged.snps_all.dbSNP142.vcf.gz
### replace chromosome name with chr ###
zcat inputs/vcf/mgp.v5.merged.snps_all.dbSNP142.vcf.gz | head -1000 | grep "^#" | cut -f 1-8 | sed 's/ID=\([0-9]\+\),/ID=chr\1,/' | sed 's/ID=\([XY]\),/ID=chr\1,/' | sed 's/ID=MT,/ID=chrM,/' > inputs/vcf/mgp.v5.merged.snps_all.dbSNP142.pass.chr.vcf
zcat inputs/vcf/mgp.v5.merged.snps_all.dbSNP142.vcf.gz | grep -v "^#" | cut -f 1-8 | grep -w "PASS" | sed 's/^/chr/' | sed 's/chrMT/chrM/' >> inputs/vcf/mgp.v5.merged.snps_all.dbSNP142.pass.chr.vcf
gatk UpdateVcfSequenceDictionary -SD inputs/fasta/mm10_AAV.dict -I inputs/vcf/mgp.v5.merged.snps_all.dbSNP142.pass.chr.vcf -O inputs/vcf/mgp.v5.merged.snps_all.dbSNP142.pass.chr.updated.vcf
gatk SortVcf -SD inputs/fasta/mm10_AAV.dict -I inputs/vcf/mgp.v5.merged.snps_all.dbSNP142.pass.chr.updated.vcf -O inputs/vcf/mgp.v5.merged.snps_all.dbSNP142.pass.chr.sort.vcf
bgzip inputs/vcf/mgp.v5.merged.snps_all.dbSNP142.pass.chr.sort.vcf
tabix inputs/vcf/mgp.v5.merged.snps_all.dbSNP142.pass.chr.sort.vcf.gz
rm inputs/vcf/mgp.v5.merged.snps_all.dbSNP142.pass.chr.vcf inputs/vcf/mgp.v5.merged.snps_all.dbSNP142.pass.chr.updated.vcf

### get known INDELs ###
wget ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz -O inputs/vcf/mgp.v5.indels.vcf.gz
### the correction of files were adapted from the post https://zqfang.github.io/2020-03-10-gatk4-mm10-bundle/ ###
### take header first ###
zcat inputs/vcf/mgp.v5.indels.vcf.gz | head -1000 | grep "^#" | cut -f 1-8 | sed 's/ID=\([0-9]\+\),/ID=chr\1,/' | sed 's/ID=\([XY]\),/ID=chr\1,/' | sed 's/ID=MT,/ID=chrM,/' > inputs/vcf/mgp.v5.indels.pass.chr.vcf
### keep only passing and append  ###
zcat inputs/vcf/mgp.v5.indels.vcf.gz | grep -v "^#" | cut -f 1-8 | grep -w "PASS" | sed 's/^/chr/' | sed 's/chrMT/chrM/' >> inputs/vcf/mgp.v5.indels.pass.chr.vcf
gatk UpdateVcfSequenceDictionary -SD inputs/fasta/mm10_AAV.dict -I inputs/vcf/mgp.v5.indels.pass.chr.vcf -O inputs/vcf/mgp.v5.indels.pass.chr.updated.vcf
gatk SortVcf -SD inputs/fasta/mm10_AAV.dict -I inputs/vcf/mgp.v5.indels.pass.chr.updated.vcf -O inputs/vcf/mgp.v5.indels.pass.chr.sort.vcf
bgzip inputs/vcf/mgp.v5.indels.pass.chr.sort.vcf
tabix inputs/vcf/mgp.v5.indels.pass.chr.sort.vcf.gz
rm inputs/vcf/mgp.v5.indels.pass.chr.vcf inputs/vcf/mgp.v5.indels.pass.chr.updated.vcf



