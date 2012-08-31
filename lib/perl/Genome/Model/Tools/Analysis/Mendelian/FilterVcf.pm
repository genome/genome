
package Genome::Model::Tools::Analysis::Mendelian::FilterVcf;     # rename this when you give the module file a different name <--

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

my $num_affected = my $affecteds_missing = my $unaffecteds_variant = my $affecteds_variant = my $affecteds_ambiguous = 0;
my %genotypes = ();

class Genome::Model::Tools::Analysis::Mendelian::FilterVcf {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		vcf_file	=> { is => 'Text', doc => "Input in VCF format", is_optional => 0, is_input => 1},
		output_file	=> { is => 'Text', doc => "Output file for pass-filter sites", is_optional => 1, is_input => 1},
		filtered_file	=> { is => 'Text', doc => "Output file for fail-filter sites", is_optional => 1, is_input => 1},
		control_samples	=> { is => 'Text', doc => "Comma-separated list of control sample names", is_optional => 1, is_input => 1},
		inheritance_model	=> { is => 'Text', doc => "Presumed model of mendelian inheritance [autosomal-dominant]", is_optional => 1, is_input => 1, default => 'autosomal-dominant'},
		min_coverage_to_refute	=> { is => 'Text', doc => "Minimum coverage to refute a possible variant in an affected", is_optional => 1, is_input => 1, default => 10},
		max_frequency_to_refute	=> { is => 'Text', doc => "Maximum observed variant allele frequency to refute a possible variant in an affected [5]", is_optional => 1, is_input => 1, default => 10},
		min_affected_variant	=> { is => 'Text', doc => "Minimum number of affecteds with variant to include", is_optional => 1, is_input => 1, default => 1},
		max_unaffected_variant	=> { is => 'Text', doc => "Maximum number of affecteds with variant to include", is_optional => 1, is_input => 1, default => 0},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Filters a VCF for variants that pass Mendelian rules of inheritance"                 
}

sub help_synopsis {
    return <<EOS
This command filters a VCF for variants that pass Mendelian rules of inheritance
EXAMPLE:	gmt analysis mendelian filter-vcf --vcf-file myVCF.vcf --output-file myVCF.pass.vcf --filtered-file myVCF.fail.vcf
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
	my $vcf_file = $self->vcf_file;

	my %control_sample = ();
	if($self->control_samples)
	{
		my @samples = split(/\,/, $self->control_samples);
		foreach my $sample (@samples)
		{
			$control_sample{$sample} = 1;
		}
	}
	
	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";		
	}
	
	if($self->filtered_file)
	{
		open(FILTERFILE, ">" . $self->filtered_file) or die "Can't open filtered outfile: $!\n";				
	}

	my %stats = ();
	
	
	## Parse the VCF file ##
	
	my $input = new FileHandle ($vcf_file);
	my $lineCounter = 0;
	
	my @sample_names = ();
	my $num_samples = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
#		warn "$lineCounter lines parsed...\n" if(!($lineCounter % 20000));
		
		my @lineContents = split(/\t/, $line);
		my $numContents = @lineContents;
		
		if(substr($line, 0, 1) eq '#')
		{
			print OUTFILE "$line\n" if($self->output_file);
			print FILTERFILE "$line\n" if($self->filtered_file);
			
			if(substr($line, 0, 6) eq '#CHROM')
			{
				for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
				{
					$sample_names[$colCounter] = $lineContents[$colCounter];
					$num_samples++;
				}

				$stats{'num_samples'} = $num_samples;
			}
		}
		else
		{
			my ($chrom, $position, $id, $ref, $var, $score, $filter, $info, $format) = split(/\t/, $line);
			$stats{'vcf_lines'}++;
			
			## Only take sites passing filters ##
			if($filter eq '.' || uc($filter) eq 'PASS')
			{
				$stats{'vcf_lines_pass'}++;
				
				## Only take sites with a variant allele #
				if($ref && $var && $var ne '.')
				{
					my $variant_type = "snp";
					if(length($ref) > 1 || (length($var) > 1 && !($var =~ ',')))
					{
						$variant_type = "indel";
					}
					$stats{'vcf_lines_pass_' . $variant_type}++;
					
					my $mendel_status = "PASS";
					my $num_affected_variant = 0;
					
					## Get the anticipated format for genotypes ##
					
					my @formatContents = split(/\:/, $format);
					my %genotype_column = ();
					my $numFormatContents = @formatContents;
					for(my $colCounter = 0; $colCounter < $numFormatContents; $colCounter++)
					{
						my $column_name = $formatContents[$colCounter];
						$genotype_column{$column_name} = $colCounter;
					}
					
					## Go through and get each sample genotype ##
					
					for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
					{
						my $sample_name = $sample_names[$colCounter];
						my @genotypeContents = split(/\:/, $lineContents[$colCounter]);
						my $genotype = $genotypeContents[$genotype_column{'GT'}];
						my $coverage = $genotypeContents[$genotype_column{'DP'}];
						my $filter = $genotypeContents[$genotype_column{'FT'}];
						
						## Only process the genotype if it has a value and is either unfiltered or for the control sample ##
						if($genotype ne '.' && ($filter eq 'PASS' || $filter eq '.' || $control_sample{$sample_name}))
						{
							$genotype = convert_genotype($ref, $var, $genotype) if($genotype ne '.');
	
							my $gt = "";
							if($genotype eq $ref . $ref)
							{
								$gt = "Ref";
							}
							elsif($genotype eq $ref . $var)
							{
								$gt = "Het";
							}
							elsif($genotype eq $var . $var)
							{
								$gt = "Hom";
							}
							else
							{
								$gt = "Missing";
							}
						
							if($gt ne "Missing")
							{
								if($control_sample{$sample_name})
								{
									## Control sample ##
									
									if($gt eq 'Het' || $gt eq 'Hom')
									{
										$mendel_status = "Control_Was_Variant";
									}
								}
								else
								{
									## Affected sample ##
									
									if($gt eq 'Ref' && $coverage >= $self->min_coverage_to_refute)
									{
										$mendel_status = "Affected_Was_Wildtype";
									}
									elsif($gt eq 'Het' || $gt eq 'Hom')
									{
										## Pass ##
										$num_affected_variant++;
									}
								}
							}							
						}
						else
						{
							## Genotype missing or was filtered out ##
						}

						
					}
					
					## See if we passed ##
					
					if($mendel_status eq 'PASS' && $num_affected_variant < $self->min_affected_variant)
					{
						$mendel_status = "Not_Enough_Affecteds_Variant";
					}
					
					$stats{'vcf_lines_pass_' . $variant_type . '_' . $mendel_status}++;
					
					## Print to corresponding file if pass or fail ##
					if($mendel_status eq 'PASS')
					{
						my $newline = join("\t", $chrom, $position, $id, $ref, $var, $score, $filter, $info, $format);
						for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
						{
							$newline .= "\t" . $lineContents[$colCounter];
						}
						print OUTFILE "$newline\n";
					}
					else
					{
						my $newline = join("\t", $chrom, $position, $id, $ref, $var, $score, $filter, "Fail=Mendel_$mendel_status", $format);
						for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
						{
							$newline .= "\t" . $lineContents[$colCounter];
						}

						print FILTERFILE "$newline\n";
					}
				}
			}
		}
		

	}
	
	close($input);


	foreach my $key (sort keys %stats)
	{
		print "$stats{$key}\t$key\n";
	}

	if($self->output_file)
	{
		close(OUTFILE);
	}
	
	if($self->filtered_file)
	{
		close(FILTERFILE);
	}
	
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub convert_genotype
{
	my ($ref, $var, $genotype) = @_;
	
	return("NN") if($genotype eq '.' || $genotype eq '0');
	
	return($ref . $ref) if($genotype eq '0/0');
	return($ref . $var) if($genotype eq '0/1');
	return($var . $var) if($genotype eq '1/1');
	
	if($var =~ '\,')
	{
		my @vars = split(/\,/, $var);
		
		my ($gt1, $gt2) = split(/\//, $genotype);

		if($gt1 == 0)
		{
			$genotype = $ref;
		}
		else
		{
			$genotype = $vars[$gt1 - 1];
		}

		if($gt2 == 0)
		{
			$genotype .= $ref;
		}
		else
		{
			$genotype .= $vars[$gt2 - 1];
		}
	}
	else
	{
		warn "Unable to convert $ref $var $genotype\n";		
	}

	return($genotype);
}


sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}


1;


