
package Genome::Model::Tools::Analysis::Mendelian::ReportVariants;     # rename this when you give the module file a different name <--

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
    byBamOrder
    code_to_genotype_returning_code
);

my %genotypes = ();

class Genome::Model::Tools::Analysis::Mendelian::ReportVariants {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variant_file	=> { is => 'Text', doc => "List of variants to consider (annotation format)", is_optional => 0, is_input => 1},
		affected_files	=> { is => 'Text', doc => "Consensus files for affected individuals", is_optional => 0, is_input => 1},
		unaffected_files	=> { is => 'Text', doc => "Consensus files for unaffected individuals", is_optional => 1, is_input => 1},
		inheritance_model	=> { is => 'Text', doc => "Presumed model of mendelian inheritance [autosomal-dominant]", is_optional => 1, is_input => 1, default => 'autosomal-dominant'},
		min_coverage_to_refute	=> { is => 'Text', doc => "Minimum coverage to refute a possible variant in an affected [10]", is_optional => 1, is_input => 1, default => 10},
		max_frequency_to_refute	=> { is => 'Text', doc => "Maximum observed variant allele frequency to refute a possible variant in an affected [5]", is_optional => 1, is_input => 1, default => 10},
		min_affected_variant	=> { is => 'Text', doc => "Minimum number of affecteds with variant to include", is_optional => 1, is_input => 1, default => 1},
		max_unaffected_variant	=> { is => 'Text', doc => "Maximum number of affecteds with variant to include", is_optional => 1, is_input => 1, default => 0},
		output_file	=> { is => 'Text', doc => "Output file for QC result", is_optional => 1, is_input => 1},
		print_all	=> { is => 'Text', doc => "If set to 1, prints all variants regardless of mendelian pattern", is_optional => 1, is_input => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Reports variants shared by affected individuals in a Mendelian disease pedigree"                 
}

sub help_synopsis {
    return <<EOS
This command reports variants shared by affected individuals in a Mendelian disease pedigree
EXAMPLE:	gmt analysis mendelian report-variants --annotation-file mySNPs.tier1 --affected-files sample1.cns,sample2.cns --output-file variants.report.tsv
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

	my $variant_file = $self->variant_file;

	my $affected_files = $self->affected_files;
	my $unaffected_files = $self->unaffected_files if($self->unaffected_files);	

	## Get user-defined settings ##
	my $min_affecteds_variant = $self->min_affected_variant;
	my $max_unaffecteds_variant = $self->max_unaffected_variant;
	my $min_coverage_to_refute = $self->min_coverage_to_refute;
	my $max_frequency_to_refute = $self->max_frequency_to_refute;

	my $inheritance_model = $self->inheritance_model;
	
	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";		
		open(EXCLUDEDFILE, ">" . $self->output_file . ".excluded") or die "Can't open outfile: $!\n";		
	}
	


	my %stats = ();
	
	## Build an array of affected individuals' genotypes ##
	
	my @affected_array = ();
	my $num_affected = my $num_unaffected = 0;
	
	my @affected_files = split(/\,/, $affected_files);
	my @unaffected_files = split(/\,/, $unaffected_files) if($unaffected_files);
	
	print "Loading Affected samples...\n";
	
	## Count the files of each type and print the header ##
	my $header = join("\t", "chrom", "chr_start", "chr_stop", "ref", "var", "unaffecteds_variant", "affecteds_variant", "affecteds_missing", "affecteds_ambiguous");
	foreach my $affected_file (@affected_files)
	{
		$num_affected++;
		$header .= "\t" if($header);
		$header .= "AFF:" . $affected_file . "\treads1\treads2\tfreq";
		load_consensus($affected_file);
	}

	foreach my $unaffected_file (@unaffected_files)
	{
		$num_unaffected++;
		$header .= "\t" if($header);
		$header .= "UNAFF:" . $unaffected_file . "\treads1\treads2\tfreq";
		load_consensus($unaffected_file);
	}

	print "$num_affected affected samples\n";
	print "$num_unaffected unaffected samples\n";



	## PRint the header ##
	
	print "$header\n";
	if($self->output_file)
	{
		print OUTFILE "variant\tnum_affected\tnum_unaffected\t$header\n";
		print EXCLUDEDFILE "variant\tnum_affected\tnum_unaffected\t$header\n";
	}


	## Print the variants ##


	my $input = new FileHandle ($variant_file);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		(my $chromosome, my $chr_start, my $chr_stop, my $ref, my $var) = split(/\t/, $line);

		if(length($ref) > 1 || $var eq "-")
		{
			## Undo the adjustment made when formatting deletions for annotation.
			$chr_start--;
		}

		my $variant_type = "SNP";
		
		if(length($ref) > 1 || $var eq "-")
		{
			$variant_type = "DEL";
		}
		elsif(length($var) > 1 || $ref eq "-")
		{
			$variant_type = "INS";
		}

		if($lineCounter >= 0)
		{
			$stats{'num_variants'}++;
						
			my $sample_genotype_string = "";

			## AFFECTED GENOTYPES ##
			
			my $affecteds_variant = my $affecteds_wildtype = my $affecteds_missing = my $affecteds_ambiguous = my $unaffecteds_variant = 0;
			my $affecteds_homozygous = 0;
			
			## See how many affecteds carry it ##
			
			foreach my $affected_file (@affected_files)
			{
#				my %genotypes = load_consensus($affected_file);
				my $sample_genotype = "-\t-\t-\t-";

				my $key = "$affected_file\t$chromosome\t$chr_start";
				
				if($genotypes{$key})
				{
					(my $sample_call, my $sample_reads1, my $sample_reads2, my $sample_freq) = split(/\t/, $genotypes{$key});
					my $sample_coverage = $sample_reads1 + $sample_reads2;
					my $sample_freq_numeric = $sample_freq;
					$sample_freq_numeric =~ s/\%//;

					if($sample_call ne $ref)
					{
						if($variant_type eq "SNP")
						{
							## We have a variant in this affected, so count it ##
							
							$affecteds_variant++;
							
							if($sample_call eq "A" || $sample_call eq "C" || $sample_call eq "G" || $sample_call eq "T")
							{
								$affecteds_homozygous++;
							}
						}
						else
						{
							## For indels, need indel support to call it variant ##
							
							if($sample_call =~ '/')
							{
								$affecteds_variant++;
							}
						}

					}
					else
					{
						if($sample_coverage >= $min_coverage_to_refute && $sample_freq_numeric <= $max_frequency_to_refute)
						{
							$affecteds_wildtype++;							
						}
						else
						{
							$affecteds_ambiguous++;
						}


					}

					$sample_call = code_to_genotype_returning_code($sample_call);					

					$sample_genotype = "$sample_call\t$sample_reads1\t$sample_reads2\t$sample_freq";
				}
				else
				{
					$affecteds_missing++;
				}
				$sample_genotype_string .= $sample_genotype . "\t";
			}
			


			## Check to see if it occurred in multiple affected samples ##

			if($affecteds_variant >= $min_affecteds_variant || $self->print_all)
			{
				## See how many unaffecteds carry it ##
				
				foreach my $unaffected_file (@unaffected_files)
				{
#					my %genotypes = load_consensus($unaffected_file);
					my $sample_genotype = "-\t-\t-\t-";
					
					my $key = "$unaffected_file\t$chromosome\t$chr_start";
					
					if($genotypes{$key})
					{
						(my $sample_call, my $sample_reads1, my $sample_reads2, my $sample_freq) = split(/\t/, $genotypes{$key});
	
						if($sample_call ne $ref)
						{
							if($variant_type eq "SNP")
							{
								## We have a variant in this affected, so count it ##
								
								$unaffecteds_variant++;
							}
							else
							{
								## For indels, need indel support to call it variant ##
								
								if($sample_call =~ '/')
								{
									$unaffecteds_variant++;
								}
							}

						}

						$sample_call = code_to_genotype_returning_code($sample_call);
	
						$sample_genotype = "$sample_call\t$sample_reads1\t$sample_reads2\t$sample_freq";
					}
	
					$sample_genotype_string .= $sample_genotype . "\t";
				}
	

				$stats{'multiple_affecteds'}++;

				my $include_variant_flag = 0;
				my $exclude_reason = "Unknown";

				## Proceed if we found few enough unaffecteds with the variant ##
				if($self->inheritance_model)
				{
					if($self->inheritance_model eq "autosomal-dominant")
					{
						if($unaffecteds_variant > $max_unaffecteds_variant)
						{
							$stats{'in_unaffected'}++;
							$exclude_reason = "PresentInUnaffected";
						}
						elsif($affecteds_wildtype > 0)
						{
							$stats{'affected_was_wildtype'}++;
							$exclude_reason = "AffectedWasWildtype";
						}
						elsif($affecteds_homozygous == $affecteds_variant)
						{
							$stats{'all_affecteds_homozygous'}++;
							$exclude_reason = "AffectedsHomozygous";
						}
						else
						{
							$include_variant_flag = 1;
							$stats{'included_variants'}++;

							if($affecteds_missing > 0)
							{
								$stats{'included_but_missing_affected'}++;
							}
							elsif($affecteds_ambiguous > 0)
							{
								$stats{'included_but_ambiguous_affected'}++;
							}
							else
							{
								$stats{'included_and_all_affected'}++;
							}
						}
					}
				}
				elsif($unaffecteds_variant <= $max_unaffecteds_variant || $self->print_all)
				{
					$include_variant_flag = 1;
					$stats{'included_variants'}++;

				}
				else
				{
					$exclude_reason = "Only-$unaffecteds_variant-Unaffecteds-Variant";
				}
				
				if($include_variant_flag)
				{
					print "$chromosome\t$chr_start\t$chr_stop\t$ref\t$var\t";
#					print "$affecteds_wildtype affecteds WT\t$exclude_reason\t";
					print "$unaffecteds_variant\t$affecteds_variant var, $affecteds_missing miss, $affecteds_ambiguous ambig\t";
					print "$sample_genotype_string";
					print "\n";


					if($self->output_file)
					{
						print OUTFILE "$line\t";
#						print OUTFILE "$affecteds_variant\t$unaffecteds_variant\t";
						print OUTFILE join("\t", $unaffecteds_variant, $affecteds_variant, $affecteds_missing, $affecteds_ambiguous);
						print OUTFILE "\t$sample_genotype_string";
						
						if($self->print_all)
						{
							if($unaffecteds_variant == 0 && $affecteds_variant >= 2)
							{
								print OUTFILE "\t1";
							}
						}
						
						print OUTFILE "\n";
					}					
				}
				else
				{
					if($self->output_file)
					{
						print EXCLUDEDFILE "$line\t";
						print EXCLUDEDFILE join("\t", $unaffecteds_variant, $affecteds_variant, $affecteds_missing, $affecteds_ambiguous);
						print EXCLUDEDFILE "\t$sample_genotype_string\t$exclude_reason\n";						
					}
				}
			}
			else
			{				
				print EXCLUDEDFILE "$line\t";
				print EXCLUDEDFILE join("\t", $unaffecteds_variant, $affecteds_variant, $affecteds_missing, $affecteds_ambiguous);
				print EXCLUDEDFILE "\t$sample_genotype_string\t$affecteds_variant-AffectedsVariant\n";										
			}

		}		
		
	}
	
	close($input);
	
	if($self->output_file)
	{
		close(OUTFILE);
		close(EXCLUDEDFILE);
	}
	
	print $stats{'num_variants'} . " variants\n";
	print $stats{'multiple_affecteds'} . " were present in multiple affected individuals\n";
	print $stats{'in_unaffected'} . " were present in a control individual (excluded)\n";
	print $stats{'affected_was_wildtype'} . " were wildtype in an affected individual (excluded)\n";
	print $stats{'all_affecteds_homozygous'} . " were excluded because all affecteds were homozygous\n";
	print $stats{'included_variants'} . " variants were included in output\n";
	print $stats{'included_and_all_affected'} . " were variant in all affecteds\n";
	print $stats{'included_but_missing_affected'} . " were variant or missing in all affecteds\n";
	print $stats{'included_but_ambiguous_affected'} . " were called wildtype but look ambiguous in at least one affected\n";

#	print $stats{'not_in_unaffected'} . " were NOT present in unaffected individuals\n";

	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




################################################################################################
# Load Genotypes
#
################################################################################################

sub load_consensus
{                               # replace with real execution logic.
	my $genotype_file = shift(@_);
#	my %genotypes = ();
	
	my $input = new FileHandle ($genotype_file);
	my $lineCounter = 0;
	my $gtCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		(my $chrom, my $position, my $ref, my $cns, my $reads1, my $reads2, my $var_freq) = split(/\t/, $line);

		if($ref =~ /[0-9]/)
		{
			($chrom, $position, my $stop, $ref, $cns, $reads1, $reads2, $var_freq) = split(/\t/, $line);			
		}

		if(length($ref) > 1 || length($cns) > 1 || $ref eq "-" || $cns eq "-")
		{
			## If CNS is not formatted, do so ##

			if(!($cns =~ '/'))
			{
				$cns = "$ref/$cns";				
			}

		}
#		my $key = "$chrom\t$position";
		my $key = "$genotype_file\t$chrom\t$position";
		if($cns ne "N")
		{
			$genotypes{$key} = "$cns\t$reads1\t$reads2\t$var_freq";				
		}

	}
	close($input);
                            # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)

#	return(%genotypes);
}

1;
