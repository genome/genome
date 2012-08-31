
package Genome::Model::Tools::Analysis::LaneQc::CompareLoh;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# SearchRuns - Search the database for runs
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	04/01/2009 by D.K.
#	MODIFIED:	04/01/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::LaneQc::CompareLoh {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		loh_file	=> { is => 'Text', doc => "Three-column file of genotype calls chrom, pos, genotype", is_optional => 0, is_input => 1 },
		normal_variant_file	=> { is => 'Text', doc => "Variant calls in SAMtools pileup-consensus format", is_optional => 1, is_input => 1 },
		tumor_variant_file	=> { is => 'Text', doc => "Alternatively, provide a BAM file", is_optional => 1, is_input => 1 },		
		min_depth_het	=> { is => 'Text', doc => "Minimum depth to compare a het call [8]", is_optional => 1, is_input => 1},
		output_file	=> { is => 'Text', doc => "Output file for QC result", is_optional => 1, is_input => 1},
		verbose	=> { is => 'Text', doc => "Turns on verbose output [0]", is_optional => 1, is_input => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Compare tumor and normal variant calls to LOH positions from SNP array"                 
}

sub help_synopsis {
    return <<EOS
This command compares tumor and normal variant calls to LOH positions from SNP array to identify tumor-normal switches
EXAMPLE:	gmt analysis lane-qc compare-loh --loh-file [LOH] --normal [normal.snps] --tumor [tumor.snps] --output [output]
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 

EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
	my $self = shift;

	my $loh_file = $self->loh_file;
	my $normal_file = $self->normal_variant_file;
	my $tumor_file = $self->tumor_variant_file;

	print "Loading LOH calls from $loh_file...\n" if($self->verbose);
	our %genotypes = load_genotypes($loh_file);
	
	my $min_depth_het = 8;
	$min_depth_het = $self->min_depth_het if($self->min_depth_het);
	
	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";
	}

	warn "Comparing normal variants to array LOH calls...\n";
	my $normal_result = compare_variants($normal_file, $min_depth_het);

	warn "Comparing tumor variants to array LOH calls...\n";	
	my $tumor_result = compare_variants($tumor_file, $min_depth_het);	


	my ($normal_comparisons, $normal_matched_normal, $normal_matched_tumor) = split(/\t/, $normal_result);
	my ($tumor_comparisons, $tumor_matched_normal, $tumor_matched_tumor) = split(/\t/, $tumor_result);

	my $pct_normal_matched_normal = my $pct_normal_matched_tumor = my $pct_tumor_matched_normal = my $pct_tumor_matched_tumor = "";
	
	if($normal_comparisons)
	{
		$pct_normal_matched_normal = sprintf("%.2f", ($normal_matched_normal / $normal_comparisons * 100));
		$pct_normal_matched_tumor = sprintf("%.2f", ($normal_matched_tumor / $normal_comparisons * 100));
	}

	if($tumor_comparisons)
	{
		$pct_tumor_matched_normal = sprintf("%.2f", ($tumor_matched_normal / $tumor_comparisons * 100));
		$pct_tumor_matched_tumor = sprintf("%.2f", ($tumor_matched_tumor / $tumor_comparisons * 100));
	}


	my $check_result = "Unknown";
	
	if($normal_comparisons && $tumor_comparisons)
	{
		if($pct_normal_matched_normal > $pct_tumor_matched_normal && $pct_tumor_matched_tumor > $pct_normal_matched_tumor)
		{
			$check_result = "OK";
		}
		elsif($pct_normal_matched_normal < $pct_tumor_matched_normal && $pct_tumor_matched_tumor < $pct_normal_matched_tumor)
		{
			$check_result = "SampleSwitch";
		}
	}

	print "normal_comparisons\tmatched_normal\tpct_match_normal\tmatched_tumor\tpct_match_tumor\ttumor_comparisons\tmatched_normal\tpct_match_normal\tmatched_tumor\tpct_match_tumor\tqc_result\n";
	print join("\t", $normal_comparisons, $normal_matched_normal, $pct_normal_matched_normal . '%', $normal_matched_tumor, $pct_normal_matched_tumor . '%') . "\t";
	print join("\t", $tumor_comparisons, $tumor_matched_normal, $pct_tumor_matched_normal . '%', $tumor_matched_tumor, $pct_tumor_matched_tumor . '%') . "\t";
	print "$check_result\n";

	if($self->output_file)
	{
		print OUTFILE "normal_comparisons\tmatched_normal\tpct_match_normal\tmatched_tumor\tpct_match_tumor\ttumor_comparisons\tmatched_normal\tpct_match_normal\tmatched_tumor\tpct_match_tumor\tqc_result\n";
		print OUTFILE join("\t", $normal_comparisons, $normal_matched_normal, $pct_normal_matched_normal . '%', $normal_matched_tumor, $pct_normal_matched_tumor . '%') . "\t";
		print OUTFILE join("\t", $tumor_comparisons, $tumor_matched_normal, $pct_tumor_matched_normal . '%', $tumor_matched_tumor, $pct_tumor_matched_tumor . '%') . "\t";
		print OUTFILE "$check_result\n";		
	}

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


################################################################################################
# Load Genotypes
#
################################################################################################

sub load_genotypes
{                               # replace with real execution logic.
	my $loh_file = shift(@_);
	my %genotypes = ();
	
	my $input = new FileHandle ($loh_file);
	my $lineCounter = 0;
	my $gtCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		(my $chrom, my $position, my $normal_genotype, my $tumor_genotype) = split(/\t/, $line);

		my $key = "$chrom\t$position";
		
		$genotypes{$key} = "$normal_genotype\t$tumor_genotype";
		$gtCounter++;
	}
	close($input);

#	print "$gtCounter genotypes loaded\n";
	
	return(%genotypes);                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




################################################################################################
# Execute - the main program logic
#
################################################################################################

sub compare_variants
{
	my ($variant_file, $min_depth_het) = @_;

	my %stats;
	our %genotypes;

	my $input = new FileHandle ($variant_file);
	my $lineCounter = 0;

	my $file_type = "samtools";
	my $verbose_output = "";

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my @lineContents = split(/\t/, $line);
		my $chrom = $lineContents[0];
		my $position = $lineContents[1];
		my $ref_base = $lineContents[2];
		my $cns_call = $lineContents[3];
		
		my $depth = 0;
		
		if(lc($chrom) =~ "chrom")
		{
			## Ignore header ##
			$file_type = "varscan";
		}
		else
		{
			if($lineContents[6] && $lineContents[6] =~ '%')
			{
				$file_type = "varscan";
			}

			## Only check SNP calls ##
	
			if($ref_base ne "*" && length($ref_base) == 1 && length($cns_call) == 1) #$ref_base ne $cns_call
			{
				## Get depth and consensus genotype ##
	
				my $cons_gt = "";			
	
				if($file_type eq "varscan" && $cns_call ne "A" && $cns_call ne "C" && $cns_call ne "G" && $cns_call ne "T")
				{
					## Varscan CNS format ##
					$depth = $lineContents[4] + $lineContents[5];
					$cons_gt = code_to_genotype($cns_call);			
				}
				elsif($file_type eq "varscan")
				{
					## Varscan SNP format ##
					$depth = $lineContents[4] + $lineContents[5];
					my $var_freq = $lineContents[6];
					my $allele1 = $lineContents[2];
					my $allele2 = $lineContents[3];
					$var_freq =~ s/\%//;
					if($var_freq >= 80)
					{
						$cons_gt = $allele2 . $allele2;
					}
					else
					{
						$cons_gt = $allele1 . $allele2;
						$cons_gt = sort_genotype($cons_gt);
					}					
				}
				
				else
				{
					$depth = $lineContents[7];
					$cons_gt = code_to_genotype($cns_call);
				}
	
				$stats{'num_snps'}++;
				
#				warn "$stats{'num_snps'} lines parsed...\n" if(!($stats{'num_snps'} % 10000));
	
				my $key = "$chrom\t$position";
					
				if($genotypes{$key})
				{
					$stats{'num_with_genotype'}++;
					
					if($depth >= $min_depth_het)
					{
						$stats{'num_min_depth'}++;
					
						my $match_status = "";
						
						(my $normal_genotype, my $tumor_genotype) = split(/\t/, $genotypes{$key});
						
						my $normal_gt = sort_genotype($normal_genotype);
						my $tumor_gt = sort_genotype($tumor_genotype);
						my $ref_gt = code_to_genotype($ref_base);
						
						if($cons_gt eq $normal_gt)
						{
							$stats{'num_matched_normal'}++;
							$match_status .= "MatchNormal";
						}

						if($cons_gt eq $tumor_gt)
						{
							$stats{'num_matched_tumor'}++;
							$match_status .= "MatchTumor";
						}

						$match_status = "MatchNeither" if(!$match_status);
					
						
					}
				}
			
			}

		}
		

		
	}
	
	close($input);
	
	## Set zero values ##
	
	$stats{'num_matched_normal'} = 0 if(!$stats{'num_matched_normal'});
	$stats{'num_matched_tumor'} = 0 if(!$stats{'num_matched_tumor'});

	## Calculate pct ##
	
	$stats{'pct_normal_match'} = $stats{'pct_tumor_match'} = "0.00";

	if($stats{'num_min_depth'})
	{
		$stats{'pct_normal_match'} = $stats{'num_matched_normal'} / $stats{'num_min_depth'} * 100;
		$stats{'pct_normal_match'} = sprintf("%.3f", $stats{'pct_normal_match'});

		$stats{'pct_tumor_match'} = $stats{'num_matched_tumor'} / $stats{'num_min_depth'} * 100;
		$stats{'pct_tumor_match'} = sprintf("%.3f", $stats{'pct_tumor_match'});

	}
	else
	{
		$stats{'pct_normal_match'} = $stats{'pct_tumor_match'} = "--";
	}

	return(join("\t", $stats{'num_min_depth'}, $stats{'num_matched_normal'}, $stats{'num_matched_tumor'}));
	
}



################################################################################################
# Sorting
#
################################################################################################


sub byBamOrder
{
	my ($chrom_a, $pos_a) = split(/\t/, $a);
	my ($chrom_b, $pos_b) = split(/\t/, $b);
	
	$chrom_a =~ s/X/9\.1/;
	$chrom_a =~ s/Y/9\.2/;
	$chrom_a =~ s/MT/25/;
	$chrom_a =~ s/M/25/;
	$chrom_a =~ s/NT/99/;
	$chrom_a =~ s/[^0-9\.]//g;

	$chrom_b =~ s/X/9\.1/;
	$chrom_b =~ s/Y/9\.2/;
	$chrom_b =~ s/MT/25/;
	$chrom_b =~ s/M/25/;
	$chrom_b =~ s/NT/99/;
	$chrom_b =~ s/[^0-9\.]//g;

	$chrom_a <=> $chrom_b
	or
	$pos_a <=> $pos_b;
}


################################################################################################
# Load Genotypes
#
################################################################################################

sub is_heterozygous
{
	my $gt = shift(@_);
	(my $a1, my $a2) = split(//, $gt);
	return(1) if($a1 ne $a2);
	return(0);
}



################################################################################################
# Load Genotypes
#
################################################################################################

sub is_homozygous
{
	my $gt = shift(@_);
	(my $a1, my $a2) = split(//, $gt);
	return(1) if($a1 eq $a2);
	return(0);
}



################################################################################################
# Load Genotypes
#
################################################################################################

sub flip_genotype
{
	my $gt = shift(@_);
	(my $a1, my $a2) = split(//, $gt);

	if($a1 eq "A")
	{
		$a1 = "T";
	}
	elsif($a1 eq "C")
	{
		$a1 = "G";
	}
	elsif($a1 eq "G")
	{
		$a1 = "C";
	}	
	elsif($a1 eq "T")
	{
		$a1 = "A";		
	}

	if($a2 eq "A")
	{
		$a2 = "T";
	}
	elsif($a2 eq "C")
	{
		$a2 = "G";
	}
	elsif($a2 eq "G")
	{
		$a2 = "C";
	}	
	elsif($a2 eq "T")
	{
		$a2 = "A";		
	}
	
	$gt = $a1 . $a2;
	$gt = sort_genotype($gt);
	return($gt);
}

################################################################################################
# Load Genotypes
#
################################################################################################

sub sort_genotype
{
	my $gt = shift(@_);
	(my $a1, my $a2) = split(//, $gt);

	my @unsorted = ($a1, $a2);
	my @sorted = sort @unsorted;
	$a1 = $sorted[0];
	$a2 = $sorted[1];
	return($a1 . $a2);
}



sub code_to_genotype
{
	my $code = shift(@_);
	
	return("AA") if($code eq "A");
	return("CC") if($code eq "C");
	return("GG") if($code eq "G");
	return("TT") if($code eq "T");

	return("AC") if($code eq "M");
	return("AG") if($code eq "R");
	return("AT") if($code eq "W");
	return("CG") if($code eq "S");
	return("CT") if($code eq "Y");
	return("GT") if($code eq "K");

#	warn "Unrecognized ambiguity code $code!\n";

	return("NN");	
}



sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}


1;

