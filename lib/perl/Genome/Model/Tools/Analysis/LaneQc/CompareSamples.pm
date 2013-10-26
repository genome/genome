
package Genome::Model::Tools::Analysis::LaneQc::CompareSamples;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# CompareSamples - Compare filtered SNP calls between two samples to ensure that they come from the same individual
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	11/01/2010 by D.K.
#	MODIFIED:	02/09/2011 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it
use Genome::Model::Tools::Analysis::Helpers qw(
    code_to_genotype
    commify
    is_heterozygous
    is_homozygous
    sort_genotype
);

my %stats = ();

class Genome::Model::Tools::Analysis::LaneQc::CompareSamples {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variant_file1	=> { is => 'Text', doc => "Filtered SNP calls for sample 1", is_optional => 0, is_input => 1 },
		variant_file2	=> { is => 'Text', doc => "Filtered SNP calls for sample 2", is_optional => 0, is_input => 1 },
		sample_name	=> { is => 'Text', doc => "Name for the sample to use in first output column", is_optional => 1, is_input => 1 },
		min_depth_het	=> { is => 'Text', doc => "Minimum depth to compare a het call [8]", is_optional => 1, is_input => 1},
		min_depth_hom	=> { is => 'Text', doc => "Minimum depth to compare a hom call [4]", is_optional => 1, is_input => 1},
		verbose	=> { is => 'Text', doc => "Turns on verbose output [0]", is_optional => 1, is_input => 1},
		output_file	=> { is => 'Text', doc => "Output file for QC result", is_optional => 1, is_input => 1, is_output => 1}
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Compares filtered SNP calls to between two samples"                 
}

sub help_synopsis {
    return <<EOS
This command compares filtered SNP calls between two samples and returns their concordance
EXAMPLE:	gmt analysis lane-qc compare-samples --variant-file1 sample1.filtered.snp --variant-file2 sample2.filtered.snp
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

	## Get required parameters ##
	my $sample_name = "Sample";

	my $variant_file1 = $self->variant_file1;
	my $variant_file2 = $self->variant_file2;

	$sample_name = $self->sample_name if($self->sample_name);
	my $min_depth_hom = 8;
	my $min_depth_het = 12;
	$min_depth_hom = $self->min_depth_hom if($self->min_depth_hom);
	$min_depth_het = $self->min_depth_het if($self->min_depth_het);
	
	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";
	}

	
	$stats{'num_snps'} = $stats{'num_min_depth'} = $stats{'num_with_genotype'} = $stats{'num_with_variant'} = $stats{'num_variant_match'} = 0;
	$stats{'het_was_hom'} = $stats{'hom_was_het'} = $stats{'het_was_diff_het'} = $stats{'rare_hom_match'} = $stats{'rare_hom_total'} = 0;
	$stats{'num_ref_was_ref'} = $stats{'num_ref_was_hom'} = $stats{'num_ref_was_het'} = 0;
	$stats{'num_chip_was_reference'} = 0;
	
	print "Loading variants from file 1...\n";
	my %variant_calls1 = load_variant_calls($variant_file1, $min_depth_het, $min_depth_hom);
	print $stats{'num_snps'} . " SNPs loaded\n";
	$stats{'num_snps_file1'} = $stats{'num_snps'};
	$stats{'num_snps'} = 0;
	
	print "Loading variants from file 2...\n";
	my %variant_calls2 = load_variant_calls($variant_file2, $min_depth_het, $min_depth_hom);
	print $stats{'num_snps'} . " SNPs loaded\n";
	$stats{'num_snps_file2'} = $stats{'num_snps'};
	$stats{'num_snps'} = 0;

	foreach my $key (keys %variant_calls1)
	{
		if($variant_calls2{$key})
		{
			$stats{'num_snps_compared'}++;			

			my ($ref_base, $chip_gt) = split(/\t/, $variant_calls1{$key});
			($ref_base, my $cons_gt) = split(/\t/, $variant_calls2{$key});
			
			my $ref_gt = code_to_genotype($ref_base);
			
			if($chip_gt eq $ref_gt)
			{
				$stats{'num_chip_was_reference'}++;
			
				if(uc($cons_gt) eq $ref_gt)
				{
					$stats{'num_ref_was_ref'}++;
				}
				elsif(is_heterozygous($cons_gt))
				{
					$stats{'num_ref_was_het'}++;
				}
				else
				{
					$stats{'num_ref_was_hom'}++;
				}
			}
			elsif($chip_gt ne $ref_gt)
			{
				$stats{'num_with_variant'}++;
				
				my $comparison_result = "Unknown";
				
				if(is_homozygous($chip_gt))
				{
					$stats{'rare_hom_total'}++;
				}
			
				if($chip_gt eq $cons_gt)
				{
					$stats{'num_variant_match'}++;
					if(is_homozygous($chip_gt))
					{
						$stats{'rare_hom_match'}++;
					}
					
					$comparison_result = "Match";

				}
				elsif(is_homozygous($chip_gt) && is_heterozygous($cons_gt))
				{
					$stats{'hom_was_het'}++;
					$comparison_result = "HomWasHet";
				}
				elsif(is_heterozygous($chip_gt) && is_homozygous($cons_gt))
				{
					$stats{'het_was_hom'}++;
					$comparison_result = "HetWasHom";
				}
				elsif(is_heterozygous($chip_gt) && is_heterozygous($chip_gt))
				{
					$stats{'het_was_diff_het'}++;
					$comparison_result = "HetMismatch";
				}
				
				
				
			}
		}

	}

	print $stats{'num_snps_compared'} . " SNPs called in both samples\n";

	if($stats{'num_snps_file1'} > $stats{'num_snps_file2'})
	{
		$stats{'pct_overlap'} = sprintf("%.2f", $stats{'num_snps_compared'} / $stats{'num_snps_file2'} * 100);
	}
	else
	{
		$stats{'pct_overlap'} = sprintf("%.2f", $stats{'num_snps_compared'} / $stats{'num_snps_file1'} * 100);
	}

	## Calculate pct ##
	
	$stats{'pct_overall_match'} = "0.00";
	if($stats{'num_with_variant'} || $stats{'num_chip_was_reference'})
	{
		$stats{'pct_overall_match'} = ($stats{'num_variant_match'} + $stats{'num_ref_was_ref'}) / ($stats{'num_chip_was_reference'} + $stats{'num_with_variant'}) * 100;
		$stats{'pct_overall_match'} = sprintf("%.3f", $stats{'pct_overall_match'});
	}

	$stats{'pct_variant_match'} = "0.00";
	if($stats{'num_with_variant'})
	{
		$stats{'pct_variant_match'} = $stats{'num_variant_match'} / $stats{'num_with_variant'} * 100;
		$stats{'pct_variant_match'} = sprintf("%.3f", $stats{'pct_variant_match'});
	}

	$stats{'pct_hom_match'} = "0.00";
	if($stats{'rare_hom_total'})
	{
		$stats{'pct_hom_match'} = $stats{'rare_hom_match'} / $stats{'rare_hom_total'} * 100;
		$stats{'pct_hom_match'} = sprintf("%.3f", $stats{'pct_hom_match'});
	}

	if($self->verbose)
	{
		print $stats{'num_snps_compared'} . " SNPs compared\n";
		print $stats{'num_variant_match'} . " had concordant genotypes\n";
		print $stats{'hom_was_het'} . " homozygotes from file1 were called heterozygous in file2\n";
		print $stats{'het_was_hom'} . " heterozygotes from file1 were called homozygous in file2\n";
		print $stats{'het_was_diff_het'} . " heterozygotes from file1 were different heterozygote\n";
		print $stats{'pct_variant_match'} . "% concordance at variant sites\n";
		print $stats{'pct_hom_match'} . "% concordance at rare-homozygous sites\n";
		print $stats{'pct_overall_match'} . "% overall concordance match\n";
	}
	else
	{
		print "Sample\tFile1\tFile2\tOverlap\tCompared\tMatched\tConcord\n";
		print "$sample_name\t";
		print commify($stats{'num_snps_file1'}) . "\t";
		print commify($stats{'num_snps_file2'}) . "\t";
		print $stats{'pct_overlap'} . "\t";
		print commify($stats{'num_snps_compared'}) . "\t";
		print commify($stats{'num_variant_match'}) . "\t";
#		print $stats{'hom_was_het'} . "\t";
#		print $stats{'het_was_hom'} . "\t";
#		print $stats{'het_was_diff_het'} . "\t";
#		print $stats{'pct_variant_match'} . "%\t";
#		print $stats{'pct_hom_match'} . "%\t";		
		print $stats{'pct_overall_match'} . "%\n";
	}

	if($self->output_file)
	{
		print OUTFILE "Sample\tCompared\tMatched\tConcord\n";
		print OUTFILE "$sample_name\t";
		print OUTFILE $stats{'num_snps'} . "\t";
		print OUTFILE $stats{'num_variant_match'} . "\t";
#		print OUTFILE $stats{'hom_was_het'} . "\t";
#		print OUTFILE $stats{'het_was_hom'} . "\t";
#		print OUTFILE $stats{'het_was_diff_het'} . "\t";
#		print OUTFILE $stats{'pct_variant_match'} . "%\t";
#		print OUTFILE $stats{'pct_hom_match'} . "%\t";		
		print OUTFILE $stats{'pct_overall_match'} . "%\n";		
	}

	## return with the genotype concordance between 0 and 100 ##
	return $stats{'pct_overall_match'};                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


################################################################################################
# Load Genotypes
#
################################################################################################

sub load_variant_calls
{
	my ($FileName, $min_depth_het, $min_depth_hom) = @_;

	my %calls = ();
	
	# replace with real execution logic.
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;

	my $file_type = "samtools";

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my @lineContents = split(/\t/, $line);
		my $numContents = @lineContents;
		my $chrom = $lineContents[0];
		my $position = $lineContents[1];
		my $ref_base = $lineContents[2];
		my $cns_call = $lineContents[3];
		
		my $depth = 0;

		## Identify newer BED format of SNVs ##
		
		if($lineContents[1] =~ /[0-9]/ && $lineContents[2] =~ /[0-9]/ && $lineContents[3] =~ '/')
		{
			$file_type = "bed";
			($ref_base, $cns_call) = split(/\//, $lineContents[3]);

		}
		
		if(lc($chrom) =~ "chrom")
		{
			## Ignore header ##
			$file_type = "varscan";
		}
		

		else
		{
			if($numContents == 3)
			{
				$file_type = "array";
			}
			
			if($lineContents[6] && $lineContents[6] =~ '%')
			{
				$file_type = "varscan";
			}

			## Only check SNP calls ##
	
			if($file_type eq "array" || ($ref_base ne "*" && length($ref_base) == 1 && length($cns_call) == 1)) #$ref_base ne $cns_call
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
				elsif($file_type eq "bed")
				{
#					($ref_base, $cns_call) = split(/\//, $lineContents[3]);
					$depth = $lineContents[5];
					$cons_gt = code_to_genotype($cns_call);
				}
				elsif($file_type eq "array")
				{
					$cons_gt = $lineContents[2];
					$depth = $min_depth_het;
				}
				
				else
				{
					$depth = $lineContents[7];
					$cons_gt = code_to_genotype($cns_call);
				}
	
				$stats{'num_snps'}++;

				if($depth >= $min_depth_het)
				{
					my $snp_key = join("\t", $chrom, $position);
					$calls{$snp_key} = "$ref_base\t$cons_gt";								
				}

			}

		}
		

		
	}
	
	close($input);


	return(%calls);
                     # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

1;
