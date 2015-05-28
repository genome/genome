
package Genome::Model::Tools::Analysis::LaneQc::SnpFrequencies;     # rename this when you give the module file a different name <--

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
use Genome::Model::Tools::Analysis::Helpers qw(
    code_to_genotype
    commify
    is_heterozygous
    is_homozygous
    sort_genotype
);

class Genome::Model::Tools::Analysis::LaneQc::SnpFrequencies {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		genotype_file	=> { is => 'Text', doc => "Three-column file of genotype calls chrom, pos, genotype", is_optional => 0, is_input => 1 },
		variant_file	=> { is => 'Text', doc => "Variant calls in SAMtools pileup-consensus format", is_optional => 0, is_input => 1 },
		min_depth_het	=> { is => 'Text', doc => "Minimum depth to compare a het call [4]", is_optional => 1, is_input => 1},
		min_depth_hom	=> { is => 'Text', doc => "Minimum depth to compare a hom call [4]", is_optional => 1, is_input => 1},
		verbose	=> { is => 'Text', doc => "Turns on verbose output [0]", is_optional => 1, is_input => 1},
		output_file	=> { is => 'Text', doc => "Output file for QC result", is_optional => 1, is_input => 1}
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "In-development for SNP frequency report"                 
}

sub help_synopsis {
    return <<EOS
This command is in-development for SNP frequency report
EXAMPLE:	gt analysis lane-qc compare-snps --genotype-file affy.genotypes --variant-file lane1.var
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
	my $genotype_file = $self->genotype_file;
	my $variant_file = $self->variant_file;
	my $min_depth_hom = 4;
	my $min_depth_het = 8;
	$min_depth_hom = $self->min_depth_hom if($self->min_depth_hom);
	$min_depth_het = $self->min_depth_het if($self->min_depth_het);
	
	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";
		print OUTFILE "chrom\tposition\tref\tarray_gt\thet_hom\tseq_gt\tvar_freq\n";

		open(OUTHET, ">" . $self->output_file . ".Het") or die "Can't open outfile: $!\n";
		print OUTHET "chrom\tposition\tref\tarray_gt\thet_hom\tseq_gt\tvar_freq\n";

		open(OUTHOM, ">" . $self->output_file . ".Hom") or die "Can't open outfile: $!\n";
		print OUTHOM "chrom\tposition\tref\tarray_gt\thet_hom\tseq_gt\tvar_freq\n";
	}

	
	my %stats = ();
	$stats{'num_snps'} = $stats{'num_min_depth'} = $stats{'num_with_genotype'} = $stats{'num_with_variant'} = $stats{'num_variant_match'} = 0;
	$stats{'het_was_hom'} = $stats{'hom_was_het'} = $stats{'het_was_diff_het'} = $stats{'rare_hom_match'} = $stats{'rare_hom_total'} = 0;
	$stats{'num_ref_was_ref'} = $stats{'num_ref_was_hom'} = $stats{'num_ref_was_het'} = 0;

	print "Loading genotypes from $genotype_file...\n" if($self->verbose);
	my %genotypes = load_genotypes($genotype_file);

	print "Parsing variant calls in $variant_file...\n" if($self->verbose);

	my $input = new FileHandle ($variant_file);
	my $lineCounter = 0;

	my $file_type = "samtools";



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

			## Get depth and consensus genotype ##

			my $cons_gt = "";			
			my $var_freq = "";

			if($file_type eq "varscan")
			{
				$depth = $lineContents[4] + $lineContents[5];
				$var_freq = $lineContents[6];
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
	
			## Only check SNP calls ##
	
			if($ref_base ne $cns_call && $ref_base ne "*" && length($ref_base) == 1 && length($cns_call) == 1)
			{
				$stats{'num_snps'}++;
	
				my $key = "$chrom\t$position";
				
				if($genotypes{$key})
				{
					$stats{'num_with_genotype'}++;
					
					my $chip_gt = sort_genotype($genotypes{$key});

					if((is_homozygous($chip_gt) && $depth >= $min_depth_hom) || (is_heterozygous($chip_gt) && $depth >= $min_depth_het))
					{
						my $ref_gt = code_to_genotype($ref_base);
												
						$stats{'num_min_depth'}++;
					
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
								## Homozygous according to chip data ##
								$stats{'rare_hom_total'}++;
								if($self->output_file)
								{
									print OUTFILE "$chrom\t$position\t$ref_base\t$chip_gt\tHom\t$cons_gt\t$var_freq\n";
									print OUTHOM "$chrom\t$position\t$ref_base\t$chip_gt\tHom\t$cons_gt\t$var_freq\n";
								}
							}
							else
							{
								## Heterozygous according to chip data ##
								if($self->output_file)
								{
									print OUTFILE "$chrom\t$position\t$ref_base\t$chip_gt\tHet\t$cons_gt\t$var_freq\n";
									print OUTHET "$chrom\t$position\t$ref_base\t$chip_gt\tHet\t$cons_gt\t$var_freq\n";
								}
							}
						
							if($chip_gt eq $cons_gt)
							{
								$stats{'num_variant_match'}++;
								if(is_homozygous($chip_gt))
								{
									$stats{'rare_hom_match'}++;
								}
								else
								{
									## Match at het site ##
									if($var_freq)
									{
#										print "$chrom\t$position\t$ref_base\t$chip_gt\t$cons_gt\t$var_freq\n";
									}									
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
							
							if($self->verbose)
							{
								print "$line\t$chip_gt $comparison_result $cons_gt\n";
							}
							
							
						}
					}
				}
			
			}			
		}
		

		
	}
	
	close($input);

	## Parse out info from variant file ##

	my @fileContents = split(/\//, $variant_file);
	my $numContents = @fileContents;
	my $lane_info = $fileContents[$numContents - 2];
	my $machine_info = $fileContents[$numContents - 3];
	my @machineContents = split(/\_/, $machine_info);
	$numContents = @machineContents;
	my $flowcell = $machineContents[$numContents - 1];
	(my $lane) = split(/\_/, $lane_info);


	## Calculate pct ##
	
	$stats{'pct_overall_match'} = "0.00";
	if($stats{'num_with_variant'} || $stats{'num_chip_was_reference'})
	{
		$stats{'pct_overall_match'} = ($stats{'num_variant_match'}) / ($stats{'num_chip_was_reference'} + $stats{'num_with_variant'}) * 100;
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
		print $stats{'num_snps'} . " SNPs parsed from variants file\n";
		print $stats{'num_with_genotype'} . " had genotype calls from the SNP array\n";
		print $stats{'num_min_depth'} . " met minimum depth of >= $min_depth_hom/$min_depth_het\n";
		print $stats{'num_chip_was_reference'} . " were called Reference on chip\n";
#		print $stats{'num_ref_was_ref'} . " reference were called reference\n";
		print $stats{'num_ref_was_het'} . " reference were called heterozygous\n";
		print $stats{'num_ref_was_hom'} . " reference were called homozygous\n";
		print $stats{'num_with_variant'} . " had informative genotype calls\n";
		print $stats{'num_variant_match'} . " had matching calls from sequencing\n";
		print $stats{'hom_was_het'} . " homozygotes from array were called heterozygous\n";
		print $stats{'het_was_hom'} . " heterozygotes from array were called homozygous\n";
		print $stats{'het_was_diff_het'} . " heterozygotes from array were different heterozygote\n";
		print $stats{'pct_variant_match'} . "% concordance at variant sites\n";
		print $stats{'pct_hom_match'} . "% concordance at rare-homozygous sites\n";
		print $stats{'pct_overall_match'} . "% overall concordance match\n";
	}
	else
	{
		print "$flowcell.$lane\t";
		print $stats{'num_snps'} . "\t";
		print $stats{'num_with_genotype'} . "\t";
		print $stats{'num_min_depth'} . "\t";
		print $stats{'num_chip_was_reference'} . "\t";
#		print $stats{'num_ref_was_ref'} . "\t";
		print $stats{'num_ref_was_het'} . "\t";
		print $stats{'num_ref_was_hom'} . "\t";
		print $stats{'num_with_variant'} . "\t";
		print $stats{'num_variant_match'} . "\t";
		print $stats{'hom_was_het'} . "\t";
		print $stats{'het_was_hom'} . "\t";
		print $stats{'het_was_diff_het'} . "\t";
		print $stats{'pct_variant_match'} . "%\t";
		print $stats{'pct_hom_match'} . "%\t";		
		print $stats{'pct_overall_match'} . "%\n";
	}

	if($self->output_file)
	{
		print OUTFILE "$flowcell.$lane\t";
		print OUTFILE $stats{'num_snps'} . "\t";
		print OUTFILE $stats{'num_with_genotype'} . "\t";
		print OUTFILE $stats{'num_min_depth'} . "\t";
		print OUTFILE $stats{'num_chip_was_reference'} . "\t";
#		print OUTFILE $stats{'num_ref_was_ref'} . "\t";
		print OUTFILE $stats{'num_ref_was_het'} . "\t";
		print OUTFILE $stats{'num_ref_was_hom'} . "\t";
		print OUTFILE $stats{'num_with_variant'} . "\t";
		print OUTFILE $stats{'num_variant_match'} . "\t";
		print OUTFILE $stats{'hom_was_het'} . "\t";
		print OUTFILE $stats{'het_was_hom'} . "\t";
		print OUTFILE $stats{'het_was_diff_het'} . "\t";
		print OUTFILE $stats{'pct_variant_match'} . "%\t";
		print OUTFILE $stats{'pct_hom_match'} . "%\t";		
		print OUTFILE $stats{'pct_overall_match'} . "%\n";		
	}

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


################################################################################################
# Load Genotypes
#
################################################################################################

sub load_genotypes
{                               # replace with real execution logic.
	my $genotype_file = shift(@_);
	my %genotypes = ();
	
	my $input = new FileHandle ($genotype_file);
	my $lineCounter = 0;
	my $gtCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		(my $chrom, my $position, my $genotype) = split(/\t/, $line);

		my $key = "$chrom\t$position";
		
		if($genotype && $genotype ne "--")
		{
			$genotypes{$key} = $genotype;
			$gtCounter++;
		}
	}
	close($input);

#	print "$gtCounter genotypes loaded\n";
	
	return(%genotypes);                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

1;
