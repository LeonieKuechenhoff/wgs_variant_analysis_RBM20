#### Author: Shengdi Li ####
#### 2022 July. 18 ####
#### Extracting attribute tables from a vcf file ####
#### The script is designed for extracting sub-columns in a vcf files and transform into tab-delimited tables. ####

use Getopt::Long;
use Bio::SeqIO;
use strict;

GetOptions('h' => \my $display_help, 'in=s' => \my $invcf, 'out=s' =>\my $outtable, 's=s' =>\my @attributes);
if($display_help)
{
	print "Command: \n\tperl vcf_to_AFtable.pl [-h] -in INPUT_VCF -out OUTPUT_TABLE\n\n";
	print "Function: \n\tLoad input file in vcf format and output allele frequency table...\n\n";
	print "Usage: \n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tInput list file containing the meta data of samples.\n";
	print "\t-out\tOutput directory to put all result subfolders.\n";
	print "\t-s\tNames of attributes and samples to be outputed (format ex.: AD|H279,L279,T279).\n";
	exit 1;
}

if((!$invcf)||(!$outtable)||(!@attributes))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}

#### load sample names into hash ####
my %samples;
foreach my $y (@attributes)
{
	my ($filter,$sample_names)=split("\\|",$y);
	foreach my $x (split(",",$sample_names))
	{
		$samples{$filter}{$x}++;
	}
}

### load vcf file ####
### check if is vcf.gz file ###
if($invcf =~ /vcf\.gz$/)
{
	open(IN,"zcat ".$invcf." | ");
}
else
{
	open(IN,$invcf);
}
open(OUT,">".$outtable);
print OUT "chr\tpos\tref\talt";
foreach my $y (@attributes)
{
	my ($filter,$sample_names)=split("\\|",$y);
	foreach my $x (split(",",$sample_names))
	{
		print OUT "\t",$filter."|".$x;
	}
}
print OUT "\n";
### read per line ###
while(my $read_line=<IN>)
{
	chomp $read_line;
	### if is annotation line ###
	if($read_line=~/^#/)
	{
		if($read_line=~/^#CHROM/)
		{
			my @read_column=split("\t",$read_line);
			my $len=@read_column;
			for(my $i=0;$i<$len;$i++)
			{
				foreach my $y (@attributes)
				{
					my ($filter,$sample_names)=split("\\|",$y);
					if($samples{$filter}{$read_column[$i]})
					{
						$samples{$filter}{$read_column[$i]} = $i
					}
				}
			}
		}
	}
	### if is variant line ###
	else
	{
		my @read_column=split("\t",$read_line);
		#### print chromosome, position, ref, alt to file ####
		print OUT $read_column[0],"\t",$read_column[1],"\t",$read_column[3],"\t",$read_column[4];
		my %column_index;
		my @split_INFO=split(":",$read_column[8]);
		my $tmplen=@split_INFO;
		for(my $j=0;$j<$tmplen;$j++)
		{
			foreach my $y (@attributes)
			{
				my ($filter,$sample_names)=split("\\|",$y);
				if($split_INFO[$j] eq $filter)
				{
					$column_index{$filter} = $j;
				}
			}
		}
		foreach my $y (@attributes)
		{
			my ($filter,$sample_names)=split("\\|",$y);
			foreach my $x (split(",",$sample_names))
			{
				### split the sample information into sub columns ###
				my @split_info = split(":",$read_column[$samples{$filter}{$x}]);
				print OUT "\t",$split_info[$column_index{$filter}];
			}
			
		}
		print OUT "\n";
	}
}

