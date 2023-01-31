#### Author: Shengdi Li ####
#### 2022 Aug. 29 ####
#### Performing semi-global alignment between two sequences (guide RNA vs. genome DNA) ####

use Getopt::Long;
use Bio::SeqIO;
use strict;

GetOptions('h' => \my $display_help, 'in=s' => \my $input_table, 'out=s' =>\my $output_table);
if($display_help)
{
	print "Command: \n\tperl semi_global_alignment.pl [-h] -in INPUT_TABLE -out OUTPUT_TABLE\n\n";
	print "Function: \n\tLoad input sequence pairs to be aligned. Format: seq1:gRNA+PAM; seq2:genome DNA...\n\n";
	print "Usage: \n";
	print "\t-h\tPrint help info.\n";
	print "\t-in\tInput list of gRNA+PAM and genome DNA pairs.\n";
	print "\t-out\tOutput path.\n";
	exit 1;
}



if((!$input_table)||(!$output_table))
{
	print "Missing arguments.\n";
	print "Use -h for looking at usage...\n";
	exit 1;
}

### Define nucleotide match table ###
my $match = 1;
my $mismatch = -3;
my $gop = -5;
my $gep = -10000;
my %match_score;
my @nucleotides = qw{A T C G};
foreach my $x(@nucleotides)
{
	foreach my $y(@nucleotides)
	{
		if($x eq $y)
		{
			$match_score{$x}{$y} = $match;
		}
	}
}
### Ambiguous bases ###
foreach my $x(@nucleotides)
{
	$match_score{"N"}{$x} = $match;
	$match_score{$x}{"N"} = $match;	
}
foreach my $x(qw{A T})
{
	$match_score{"W"}{$x} = $match;
	$match_score{$x}{"W"} = $match;
}
foreach my $x(qw{C G})
{
	$match_score{"S"}{$x} = $match;
	$match_score{$x}{"S"} = $match;
}
foreach my $x(qw{A C})
{
	$match_score{"M"}{$x} = $match;
	$match_score{$x}{"M"} = $match;	
}
foreach my $x(qw{G T})
{
	$match_score{"K"}{$x} = $match;
	$match_score{$x}{"K"} = $match;	
}
foreach my $x(qw{A G})
{
	$match_score{"R"}{$x} = $match;
	$match_score{$x}{"R"} = $match;	
}
foreach my $x(qw{C T})
{
	$match_score{"Y"}{$x} = $match;
	$match_score{$x}{"Y"} = $match;
}
foreach my $x(qw{C G T})
{
	$match_score{"B"}{$x} = $match;
	$match_score{$x}{"B"} = $match;
}
foreach my $x(qw{A G T})
{
	$match_score{"D"}{$x} = $match;
	$match_score{$x}{"D"} = $match;
}
foreach my $x(qw{A C T})
{
	$match_score{"H"}{$x} = $match;
	$match_score{$x}{"H"} = $match;
}
foreach my $x(qw{A C G})
{
	$match_score{"V"}{$x} = $match;
	$match_score{$x}{"V"} = $match;
}

### Define semi-global alignment ###
sub SG_align()
{
	my ($seq1,$seq2)=@_;
	### Define score matrix and trace-back matrix
	my @smat;
	my @tb;
	my $len1=length($seq1);
	my $len2=length($seq2);
	
	$tb[0][0]=-100;
	$smat[0][0]=0;
	for (my $i=1; $i<=$len1; $i++)
	{
		$smat[$i][0]=-10000;
		$tb[$i][0]= 1;
	}
	for (my $j=1; $j<=$len2; $j++)
	{
		$smat[0][$j]=$smat[0][$j-1];
		$tb[0][$j]=-1;
	}

	### Fill the score matrix ###
	### first part of sequence 1 ###
	for (my $i=1; $i<=($len1-3); $i++)
	{
		for (my $j=1; $j<=$len2; $j++)
		{
			### match/mismatch ###
			my $base1=substr($seq1,$i-1,1);
			my $base2=substr($seq2,$j-1,1);
			my $sub=$smat[$i-1][$j-1]+$match;
			if(!$match_score{$base1}{$base2})
			{
				$sub=$smat[$i-1][$j-1]+$mismatch;
			}
			### deletion ###
			my $del=$smat[$i][$j-1]+$gep;
			if($tb[$i][$j-1]!=-1)
			{
				$del=$smat[$i][$j-1]+$gop+$gep;
			}
			### insertion ###
			my $ins=$smat[$i-1][$j]+$gep;;
			if($tb[$i-1][$j]!=1)
			{
				$ins=$smat[$i-1][$j]+$gop+$gep;
			}
			
			### max(substitution, deletion, insertion)
			if(($sub>$del) && ($sub>$ins))
			{
				$smat[$i][$j]=$sub;
				$tb[$i][$j]=0;
			}
			elsif($del>$ins)
			{
				$smat[$i][$j]=$del;
				$tb[$i][$j]=-1;
			}
			else 
			{
				$smat[$i][$j]=$ins;
				$tb[$i][$j]=1;
			}
		}
	}

	### last 3 bases of seq1 (PAM) ###
	for (my $i=($len1-2); $i<=$len1; $i++)
	{
		for (my $j=1; $j<=$len2; $j++)
		{
			### match/mismatch ###
			my $base1=substr($seq1,$i-1,1);
			my $base2=substr($seq2,$j-1,1);
			my $sub=$smat[$i-1][$j-1]+$match;
			if(!$match_score{$base1}{$base2})
			{
				$sub=-10000;
			}
			### deletion in seq 2 ###
			my $del=-10000;
			### semi-global alignment allows gaps at edge of seq2 ###
			if($i==$len1)
			{
				$del=$smat[$i][$j-1];
			}
			### deletion in seq 1 ###
			### deletion in PAM sites not allowed ###
			my $ins=-10000;
				
			### max(substitution, deletion, insertion)
			if(($sub>$del) && ($sub>$ins))
			{
				$smat[$i][$j]=$sub;
				$tb[$i][$j]=0;
			}
			elsif($del>$ins)
			{
				$smat[$i][$j]=$del;
				$tb[$i][$j]=-1;
			}
			else 
			{
				$smat[$i][$j]=$ins;
				$tb[$i][$j]=1;
			}
		}
	}

	### trace back ###
	my $i=$len1;
	my $j=$len2;
	my $final_score = $smat[$i][$j];
	my $aln_len=0;
	my $mm=0;
	my @aln1;
	my @aln2;

	while (!($i==0 && $j==0))
	{
		if ($tb[$i][$j]==0)
		{
			$aln1[$aln_len]=substr($seq1,$i-1,1);
			$aln2[$aln_len]=substr($seq2,$j-1,1);
			$i--;
			$j--;
		}
		elsif ($tb[$i][$j]==-1)
		{
			$aln1[$aln_len]='-';
			$aln2[$aln_len]=substr($seq2,$j-1,1);
			$j--;
			
			
		}
		elsif ($tb[$i][$j]==1)
		{
			$aln1[$aln_len]=substr($seq1,$i-1,1);
			$aln2[$aln_len]='-';
			$i--;
		}
		$aln_len++;
		
	}
	
	my $aln1=join("",@aln1),"\n";
	my $aln2=join("",@aln2),"\n";
	$aln1=reverse $aln1;
	$aln2=reverse $aln2;
	### output the fasta ###
	my $mm = 0;
	my $del1 = 0;
	my $del2 = 0;
	my $tmp1 = $aln1;
	my $tmp2 = $aln2;
	my $matched_seq;

	while(substr($tmp1,0,1) eq "-")
	{
		substr($tmp1,0,1)="";
		substr($tmp2,0,1)="";
	}

	
	while(substr($tmp1,length($tmp1)-1,1) eq "-")
	{
		substr($tmp1,length($tmp1)-1,1)="";
		substr($tmp2,length($tmp1)-1,1)="";
	}

	
	for(my $n=0;$n<length($tmp1);$n++)
	{
		if(substr($tmp1,$n,1) eq "-")
		{
			$matched_seq=$matched_seq."{".lc(substr($tmp2,$n,1))."}";
			$del1++;
		}
		elsif(substr($tmp2,$n,1) eq "-")
		{
			$matched_seq=$matched_seq."-";
			$del2++;
		}
		elsif($match_score{substr($tmp1,$n,1)}{substr($tmp2,$n,1)})
		{
			$matched_seq= $matched_seq.substr($tmp1,$n,1);
		}
		else
		{
			$matched_seq= $matched_seq.lc(substr($tmp2,$n,1));
			$mm++;
		}
	}
	return($aln1,$aln2,$matched_seq,$final_score,$mm,$del1,$del2);
}

sub reverseComp()
{
	my $seq = $_[0];
	my %rc=("A"=>"T","T"=>"A","C"=>"G","G"=>"C",
			"N"=>"N","W"=>"S","S"=>"W","M"=>"K","K"=>"M","R"=>"Y","Y"=>"R",
			"B"=>"V","V"=>"B","D"=>"H","H"=>"D","-"=>"-");
	$seq = reverse $seq;
	for(my $i=0;$i<length($seq);$i++)
	{
		substr($seq,$i,1)=$rc{substr($seq,$i,1)}
	}
	return $seq;
}

open(IN,$input_table);
open(OUT,">".$output_table);
my $header = <IN>;
chomp $header;
print OUT $header."\tmatched_strand\tmatched_sequence\talignment_score\taln1\taln2\tmismtaches\tseq1_deletion\tseq2_deletion\n";
while(my $read_line=<IN>)
{
	chomp $read_line;
	my ($seq1,$seq2)=split("\t",$read_line);
	my $strand = "+";
	### forward strand ###
	my ($aln1,$aln2,$matched_seq,$final_score,$mm,$del1,$del2)=&SG_align($seq1,$seq2);
	my ($xaln1,$xaln2,$xmatched_seq,$xfinal_score,$xmm,$xdel1,$xdel2)=&SG_align($seq1,&reverseComp($seq2));
	if($final_score <= $xfinal_score)
	{
		$strand = "-";
		($aln1,$aln2,$matched_seq,$final_score,$mm,$del1,$del2)=&SG_align($seq1,&reverseComp($seq2));
	}
	print OUT $read_line,"\t",$strand,"\t",$matched_seq,"\t",$final_score,"\t",$aln1,"\t",$aln2,"\t",$mm,"\t",$del1,"\t",$del2,"\n";
}


