package Genome::Model::Tools::Vcf::ConvertToBeagle;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# MutationRate - Calculate the mutation rate (per megabase) given a list of mutations (e.g. tier1 SNVs) and a set of regions (e.g. coding space)
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	04/22/2011 by D.K.
#	MODIFIED:	04/22/2011 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

## Declare global statistics hash ##

my %stats = ();

class Genome::Model::Tools::Vcf::ConvertToBeagle {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		vcf_file	=> { is => 'Text', doc => "Input VCF File" , is_optional => 0},
		output_file	=> { is => 'Text', doc => "Output BEAGLE file" , is_optional => 0},
		markers_file	=> { is => 'Text', doc => "Output BEAGLE markers file" , is_optional => 0},
		chromosome	=> { is => 'Text', doc => "Output results for a single chromosome" , is_optional => 1},
		complete_only	=> { is => 'Text', doc => "Only output if all samples have genotypes" , is_optional => 1},
		min_read_depth	=> { is => 'Text', doc => "Minimum VCF read depth to accept backfilled wildtype genotypes" , is_optional => 0, default => 20},		
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Converts VCF files into BEAGLE unphased genotype files"                 
}

sub help_synopsis {
    return <<EOS
This command converts VCF files into BEAGLE unphased genotype files
EXAMPLE:	gmt vcf convert-to-beagle --vcf-file my.vcf --output-file my.beagle
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
	my $output_file = $self->output_file;
	my $markers_file = $self->markers_file;

	## Open output file ##
	
	open(OUTFILE, ">$output_file") or die "Can't open output file: $!\n";
	open(MARKERS, ">$markers_file") or die "Can't open output file: $!\n";	
        
	my @sample_names = ();
	my $numSamples = 0;

	my $numVariants = my $numVariantsPassed = my $numVariantsPassedBiallelic = 0;

	## Parse the file ##

	my $input = new FileHandle ($vcf_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		if(substr($line, 0, 1) eq '#')
		{
			## HEader lines. Ignore unless samples ##
			
			if(substr($line, 0, 6) eq '#CHROM')
			{
				## Print Beagle header line ##
				
				print OUTFILE "I\tID";
				
				## Get the sample names ##
				my @lineContents = split(/\t/, $line);
				my $numContents = @lineContents;
				
				for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
				{
					my $sample = $lineContents[$colCounter];
					$sample_names[$numSamples] = $sample;
					$numSamples++;

					print OUTFILE "\t$sample\t$sample";
				}
				
				print OUTFILE "\n";
				
				print "$numSamples samples in VCF file\n";
			}
		}
		else
		{
	                my ($chrom, $position, $id, $ref, $alt, $qual, $filter, $info, $format) = split(/\t/, $line);
			my @lineContents = split(/\t/, $line);
			my $numContents = @lineContents;
			
			if(!$self->chromosome || $chrom eq $self->chromosome)
			{
				$numVariants++;
				
				if($filter eq '.' || $filter eq 'PASS')
				{
					$numVariantsPassed++;
	
					my @alt = split(/\,/, $alt);
					my $num_alts = @alt;
	
					if($num_alts == 1)
					{
						$numVariantsPassedBiallelic++;
						my $beagle_line = "M\t" . $chrom . ":" . $position;
						
						my $chrNumber = $chrom;
						$chrNumber = 23 if($chrom eq "X");
						$chrNumber = 24 if($chrom eq "Y");
						$chrNumber = 25 if($chrom eq "MT");
						$chrNumber =~ s/[^0-9]//g;
						
						if($chrNumber <= 24)
						{
							
							my @formatColumns = split(/\:/, $format);		
							
							for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
							{
								my $gt_string = $lineContents[$colCounter];
								my @gtContents = split(/\:/, $gt_string);
								my $numGtContents = @gtContents;
			
								## Parse out the relevant information ##
								my $genotype = "?";
								my $read_depth = 0;
								my $gt_filter_status = "";
								
								for(my $gtCounter = 0; $gtCounter < $numGtContents; $gtCounter++)
								{
									my $column_name = $formatColumns[$gtCounter];
									my $value = $gtContents[$gtCounter];
									
									$read_depth = $value if($column_name eq "DP" && $value ne ".");
									$genotype = $value if($column_name eq "GT" && $value ne ".");
									$gt_filter_status = $value if($column_name eq 'FT');
								}

								## Reset filtered variants ##
								
								if($gt_filter_status && $gt_filter_status ne "." && $gt_filter_status ne "PASS")
								{
									## Failed per-site genotype filter, so mark as missing ##
									$genotype = "?";
								}
								
								## Reset variants without enough coverage to say wildtype ##
								
								if($genotype eq "0/0" && $read_depth < $self->min_read_depth)
								{
									$genotype = "?";
								}
								
								
								## Print relevant genotype ##
								
								if($genotype && $genotype ne '?')
								{
									my ($a1, $a2) = split(/\//, $genotype);
									$a1 = code_to_allele($ref, $alt, $a1);
									$a2 = code_to_allele($ref, $alt, $a2);
									
									if($a1 && $a2)
									{
										$beagle_line .= "\t$a1\t$a2";							
									}
									else
									{
										$beagle_line .= "\t?\t?";							
									}
								}
								else
								{
									$beagle_line .= "\t?\t?";
								}
								
													
							}	## Go to next sample ##
							
							if($self->complete_only && $beagle_line =~ '\?')
							{
								## Skip ##
							}
							else
							{
								print MARKERS join("\t", $chrom . ":" . $position, $position, $ref, $alt) . "\n";						
								print OUTFILE "$beagle_line\n";								
							}

			#				return(1) if($lineCounter > 50);											
						}
	
					}
	
	
				}				
			}
			

			

		}

                
		

	}
	
	print "$numVariants variants in VCF file\n";
	print "$numVariantsPassed passed filters\n";
	print "$numVariantsPassedBiallelic were biallelic and included in output\n";
	
	
	close($input);

	close(OUTFILE);
	close(MARKERS);

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



################################################################################################
# Execute - the main program logic
#
################################################################################################

sub code_to_allele
{
	my ($ref, $alt, $code) = @_;
	
	my @alt = split(/\,/, $alt);
	
	## Empty ##
	return("?") if($code eq '.');
	
	## Reference ##
	return($ref) if($code eq "0");

	## Variant ##
	return($alt[$code - 1]);	
}


1;
