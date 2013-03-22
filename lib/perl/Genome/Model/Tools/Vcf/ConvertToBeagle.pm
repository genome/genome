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
my %distances = ();

class Genome::Model::Tools::Vcf::ConvertToBeagle {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		vcf_file	=> { is => 'Text', doc => "Input VCF File" , is_optional => 0},
		output_file	=> { is => 'Text', doc => "Output BEAGLE file" , is_optional => 0},
		markers_file	=> { is => 'Text', doc => "Output BEAGLE markers file" , is_optional => 0},
		chromosome	=> { is => 'Text', doc => "Output results for a single chromosome" , is_optional => 1},
		complete_only	=> { is => 'Text', doc => "Only output if all samples have genotypes" , is_optional => 1},
		min_read_depth	=> { is => 'Text', doc => "Minimum VCF read depth to accept backfilled wildtype genotypes" , is_optional => 0, default => 20},
		hapmap_file	=> { is => 'Text', doc => "Genetic map distances from HapMap, e.g. /gscmnt/sata809/info/medseq/dkoboldt/SNPseek2/ucsc/geneticMap/genetic_map_GRCh37_chr1.txt" , is_optional => 1},
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

	
	if($self->hapmap_file)
	{
		%distances = load_genetic_map($self->hapmap_file);
	}

	## Open output file ##
	
	open(OUTFILE, ">$output_file") or die "Can't open output file: $!\n";
	open(MARKERS, ">$markers_file") or die "Can't open output file: $!\n";	
        
	my @sample_names = ();
	my $numSamples = 0;

	my $numVariants = my $numVariantsPassed = my $numVariantsPassedBiallelic = my $numVariantsPassedMap = my $numVariantsPassedMapInferred = 0;

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
	
					my $variant_key = join("\t", $chrom, $position);

					## IF we have cM values but no distances for the variant key, get a distance ##
					if($self->hapmap_file && !$distances{$variant_key})
					{
						$distances{$variant_key} = infer_distance($chrom, $position);					
						$numVariantsPassedMapInferred++ if($distances{$variant_key});
					}
	
					## Only proceed if no HapMap file provided, or if we have a key for this variant ##
					if(!$self->hapmap_file || $distances{$variant_key})
					{
						$numVariantsPassedMap++;
						my $map_position = $position;
						
						if($distances{$variant_key})
						{
							my ($map_rate, $map_cm) = split(/\t/, $distances{$variant_key});
							$map_position = $map_cm;
						}
						
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
									print MARKERS join("\t", $chrom . ":" . $position, $map_position, $ref, $alt) . "\n";						
									print OUTFILE "$beagle_line\n";								
								}
	
				#				return(1) if($lineCounter > 50);											
							}
		
						}						
					}
					

	
	
				}				
			}
			

			

		}

                
		

	}
	
	print "$numVariants variants in VCF file\n";
	print "$numVariantsPassed passed filters\n";
	print "$numVariantsPassedMap had HapMap cM values ($numVariantsPassedMapInferred were inferred), or a HapMap file wasn't specified\n";
	print "$numVariantsPassedBiallelic were biallelic and included in output\n";
	
	
	close($input);

	close(OUTFILE);
	close(MARKERS);

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




#############################################################
# parse_file - parses the file
#
#############################################################

sub load_genetic_map
{
	my $FileName = shift(@_);

	my %map = ();

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		if($lineCounter > 1)
		{
			my ($chrom, $position, $rate, $cm) = split(/\t/, $line);
			$chrom =~ s/chr//;
			my $key = join("\t", $chrom, $position);
			
			$map{$key} = join("\t", $rate, $cm);
		}
	}
	
	close($input);	

	return(%map);
}


#############################################################
# parse_file - parses the file
#
#############################################################

sub infer_distance
{
	my ($chrom, $position) = @_;

	my $nearest_up = nearest_up($chrom, $position);
	my $nearest_down = nearest_down($chrom, $position);
	
	if($nearest_up && $nearest_down)
	{
		my ($up_rate, $up_cm) = split(/\t/, $nearest_up);
		my ($down_rate, $down_cm) = split(/\t/, $nearest_down);
		
		my $avg_rate = ($up_rate + $down_rate) / 2;
		my $avg_cm = ($up_cm + $down_cm) / 2;
		
		my $distance = join("\t", $avg_rate, $avg_cm);
		return($distance);		
	}

	return();

}



#############################################################
# parse_file - parses the file
#
#############################################################

sub nearest_up
{
	my ($chrom, $position) = @_;

	my $nearest_up = 0;
	
	for(my $offset = 0; $offset < 100000; $offset++)
	{
		my $up_key = join("\t", $chrom, $position + $offset);
		
		if($distances{$up_key})
		{
			return($distances{$up_key});
		}
	}
	return();
}


#############################################################
# parse_file - parses the file
#
#############################################################

sub nearest_down
{
	my ($chrom, $position) = @_;

	my $nearest_down = 0;
	
	for(my $offset = 0; $offset < 100000; $offset++)
	{
		my $down_key = join("\t", $chrom, $position - $offset);
		
		if($distances{$down_key})
		{
			return($distances{$down_key});
		}
	}
	return();
}


#############################################################
# parse_file - parses the file
#
#############################################################

sub old_infer_distance
{
	my ($chrom, $position) = @_;

	my $nearest_up = my $nearest_down = 0;
	
	for(my $offset = 0; $offset < 1000000; $offset++)
	{
		my $up_key = join("\t", $chrom, $position + $offset);
		my $down_key = join("\t", $chrom, $position - $offset);
		
		if($distances{$up_key})
		{
			$nearest_up = $distances{$up_key} if(!$nearest_up);
		}

		if($distances{$down_key})
		{
			$nearest_down = $distances{$down_key} if(!$nearest_down);
		}

	}

	if($nearest_up && $nearest_down)
	{
		my ($up_rate, $up_cm) = split(/\t/, $nearest_up);
		my ($down_rate, $down_cm) = split(/\t/, $nearest_down);
		
		my $avg_rate = ($up_rate + $down_rate) / 2;
		my $avg_cm = ($up_cm + $down_cm) / 2;
		
		my $distance = join("\t", $avg_rate, $avg_cm);
		return($distance);
	}
	
	return();

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
