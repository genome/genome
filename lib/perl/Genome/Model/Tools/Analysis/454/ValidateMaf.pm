
package Genome::Model::Tools::Analysis::454::ValidateMaf;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ValidateMaf - Align reads with SSAHA2 or other aligner
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	02/25/2009 by D.K.
#	MODIFIED:	02/25/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::454::ValidateMaf {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		maf_file	=> { is => 'Text', doc => "Complete MAF file to be updated" },
		status_file	=> { is => 'Text', doc => "Status file from Varscan with somatic_status in 13th column" },
		snp_status_file	=> { is => 'Text', doc => "Status file from Varscan with somatic_status in 13th column", is_optional => 1 },
		indel_status_file	=> { is => 'Text', doc => "Status file from Varscan with somatic_status in 13th column", is_optional => 1 },
		positions_file		=> { is => 'Text', doc => "Positions that were targeted for validation" },
		variant_type		=> { is => 'Text', doc => "Variant type, either SNP or INDEL"},
		output_file		=> { is => 'Text', doc => "Output file for updated MAF" },
		indel_padding	=> { is => 'Text', doc => "Position difference to allow to validate indel calls [15]", is_optional => 1},
		max_size_diff	=> { is => 'Text', doc => "Maximum size difference to validate indel calls [1]", is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Updates variants in MAF file with 454 validation status"                 
}

sub help_synopsis {
    return <<EOS
This command updates variants in a MAF file with 454 validation status
EXAMPLE:	gt analysis 454 validate-maf --maf-file original.maf --status-file varscan.status.snp --output-file updated.maf --positions-file targeted.tsv
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
	my $maf_file = $self->maf_file;
	
	
	my $status_file = $self->status_file;
	my $positions_file = $self->positions_file;
	my $output_file = $self->output_file;
	my $indel_padding = 15;
	$indel_padding = $self->indel_padding if($self->indel_padding);
	my $max_size_diff = 1;
	$max_size_diff = $self->max_size_diff if($self->max_size_diff);

	if(!(-e $maf_file))
	{
		die "Error: MAF file not found!\n";
	}

	if(!(-e $status_file))
	{
		die "Error: Status file not found!\n";
	}

	if(!(-e $positions_file))
	{
		die "Error: Positions file not found!\n";
	}


	my %stats = ();

	open(OUTFILE, ">$output_file") or die "Can't open output file: $!\n";

	## Load the target positions ##
	
	my %target_positions = load_target_positions($positions_file);


	## Load the 454 validation results ##
	my %results454 = load_454_results($status_file);	
#	my %results454snp = load_454_results($snp_status_file);	
#	my %results454indel = load_454_results($indel_status_file);		
	
	## Parse the MAF file ##
	
	my $input = new FileHandle ($maf_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		if($lineCounter == 1 && $line =~ "Chrom")
		{
			## Print the MAF header ##
			print OUTFILE "$line\n";
		}
		elsif($lineCounter == 1)
		{
			warn "Warning: No Header in MAF file!\n";
		}
		else
		{
			my @lineContents = split(/\t/, $line);			
			my $gene_name = $lineContents[0];
			my $chromosome = $lineContents[4];
			my $position = $lineContents[5];
			my $trv_type = $lineContents[8];
			my $var_type = $lineContents[9];
			my $ref_base = $lineContents[10];
			my $var_base = $lineContents[11];
			$var_base = $lineContents[12] if($lineContents[12] ne $ref_base);
			my $tumor_sample = $lineContents[15];
			my $current_val_status = $lineContents[24];
	
			my $position_key = "$chromosome\t$position";
			my $variant_key = "$chromosome\t$position";#\t$ref_base\t$var_base";#\t$ref_base\t$var_base";
			
			## Convert to patient ID ##
			my $patient_id = "";		
			my @tempArray = split(/\-/, $tumor_sample);
			if($tempArray[2])
			{
				$patient_id = "H_GP-" . $tempArray[2];
			}
	
			## Check for targeted position ##
			
			## SNP VALIDATION ##
			
			if($var_type eq "SNP" && $target_positions{$position_key})
			{
				$stats{'num_targeted_positions'}++;
				
				## Check for 454 variant validation ##
		
				if($results454{$variant_key})
				{
					$stats{'num_targeted_results'}++;
				
					## Parse the 454 variant results ##
				
					my @results = split(/\t/, $results454{$variant_key});
					my $ref = $results[2];
					my $var = $results[3];
					
					my $this_var_type = "SNP";
					
					if($var && length($var) > 1)
					{
						$this_var_type = "INS" if(substr($var, 0, 1) eq '+');
						$this_var_type = "DEL" if(substr($var, 0, 1) eq '-');
					}

					if($this_var_type eq $var_type)
					{
						$stats{'num_type_matched'}++;
						
						my $normal_gt = $results[7];
						my $tumor_gt = $results[11];
						my $validation_status = $results[12];
			
						if($current_val_status eq "unknown")
						{
							$stats{'changed from ' . $current_val_status . ' to ' . $validation_status}++;						
	#						print "$variant_key\t$var_type\t$current_val_status\t==>\t$validation_status\n$results454{$variant_key}\n" if($current_val_status ne $validation_status);
							
							## Parse out the normal and tumor alleles ##			
							(my $normal_allele1, my $normal_allele2) = split(/\//, $normal_gt);
							(my $tumor_allele1, my $tumor_allele2) = split(/\//, $tumor_gt);
				
							$lineContents[19] = $tumor_allele1;
							$lineContents[20] = $tumor_allele2;
							$lineContents[21] = $normal_allele1;
							$lineContents[22] = $normal_allele2;
							$lineContents[24] = $validation_status;
						}
	
						## Build a new line for MAF ##	
						my $newline = "";
						my $numContents = @lineContents;
						for(my $colCounter = 0; $colCounter < $numContents; $colCounter++)
						{
							$lineContents[$colCounter] = "" if(!$lineContents[$colCounter]);
							$newline .= $lineContents[$colCounter] . "\t";
						}
						
						print OUTFILE "$newline\n";										
					}
					else
					{
#						print "Can't validate $variant_key $var_type with $this_var_type $results454{$variant_key}\n";
						## Variant type did not match ##
						print OUTFILE "$line\n";				
					}

				
				}
				else
				{
					## Variant had no result - print original line ##
					print OUTFILE "$line\n";				
				}
			}
			
			## INDEL VALIDATION ##
			
			elsif($current_val_status eq "unknown" && $self->variant_type eq "INDEL" && ($var_type eq "INS" || $var_type eq "DEL"))
			{

				## Determine length of region ##
				
				my $indel_size = my $region_start = my $region_stop = 0;
				my $indel_bases = "";
				
				$indel_bases = $var_base;
				$indel_bases =~ s/[^ACGTN]//g;
				$indel_size = length($indel_bases);
				
#				$region_start = $position - $indel_size;
#				$region_stop = $position + $indel_size;

				$region_start = $position - $indel_size - $indel_padding;
				$region_stop = $position + $indel_size + $indel_padding;

				## Reset variables for tracking newlines ##
				my $newline = my $new_val_status = my $new_val_position = "";
				

				for(my $this_pos = $region_start; $this_pos <= $region_stop; $this_pos++)
				{
					$position_key = "$chromosome\t$this_pos";
					if($results454{$position_key})
					{
						$stats{'indels_with_results'}++;
						
						## Parse the 454 variant results ##
					
						my @results = split(/\t/, $results454{$position_key});
						
						my $ref = $results[2];
						my $var = $results[3];
						my $normal_gt = $results[7];
						my $tumor_gt = $results[11];
						my $validation_status = $results[12];						

#						$validation_bases++;
						
						if($validation_status && $validation_status eq "Reference")
						{
#							$non_indel_bases++;
						}
						else
						{
							my $indel_type_454 = my $indel_size_454 = my $indel_bases_454 = "";
							
							$indel_bases_454 = uc($var);
							$indel_bases_454 =~ s/[^ACGTN]//g;
							
							if(substr($var, 0, 1) eq "+")
							{
								$indel_type_454 = "INS";
							}
							else
							{
								$indel_type_454 = "DEL";
							}
							
							$indel_size_454 = length($indel_bases_454);
							
							## Ensure that 454 indel type matches MAF indel type ##
							
							if($indel_type_454 eq $var_type)
							{
								$stats{'num_type_matched'}++;
								
								## Determine size diff ##
								
								my $size_diff = abs($indel_size_454 - $indel_size);
	
								## Ensure that indels are roughly the same size ##
								
								if($size_diff <= $max_size_diff)
								{
									$stats{'num_size_matched'}++;
		
									## Parse out the normal and tumor alleles ##			
									(my $normal_allele1, my $normal_allele2) = split(/\//, $normal_gt);
									(my $tumor_allele1, my $tumor_allele2) = split(/\//, $tumor_gt);
						
									$lineContents[19] = $tumor_allele1;
									$lineContents[20] = $tumor_allele2;
									$lineContents[21] = $normal_allele1;
									$lineContents[22] = $normal_allele2;
									$lineContents[24] = $validation_status;
	
									## Determine if we need to update the newline. Only do so if: ##
									# 1. This is the first validation results for this variant, or
									# 2. The previous validation result called it Reference/WildType, or
									# 3. This validation result [matches the previous one, but this] is physically closer to the MAF variant
									
									if(!$new_val_status || ($validation_status eq "Somatic" && $new_val_status ne "Somatic") || $new_val_status eq "Reference" || (abs($position - $this_pos) < abs($position - $new_val_position)))  #$validation_status eq $new_val_status 
#									if(1)
									{
										print "$gene_name\t$chromosome\t$position\t$var_type\t$ref_base\t$var_base\t$tumor_sample\t$current_val_status\n";
										print "454data\t$position_key\t$ref\t$var\tN=$normal_gt\tT=$tumor_gt\t$validation_status\t***UPDATED***\n\n";
										## Update the newline ##
	
										$new_val_status = $validation_status;
										$new_val_position = $this_pos;
	
										## Build a new line for MAF ##	
										$newline = "";
										my $numContents = @lineContents;
										for(my $colCounter = 0; $colCounter < $numContents; $colCounter++)
										{
											$lineContents[$colCounter] = "" if(!$lineContents[$colCounter]);
											$newline .= $lineContents[$colCounter] . "\t";
										}
	
										## Print the update that was made ##
#										print "$position_key\t$var_type-$indel_size\t$ref_base\t$indel_bases\n";
#										print "Updated\t$this_pos\t$indel_type_454-$indel_size_454\t$ref\t$var\t$validation_status\n";															
										
										$stats{'changed from ' . $current_val_status . ' to ' . $validation_status}++;	
									}
									else
									{
#$										print "$position_key status NOT updated from $new_val_status to $validation_status\n";
										print "$gene_name\t$chromosome\t$position\t$var_type\t$ref_base\t$var_base\t$tumor_sample\t$current_val_status\n";
										print "454data\t$position_key\t$ref\t$var\tN=$normal_gt\tT=$tumor_gt\t$validation_status\t***NOT CHANGED***\n\n";
									}
	
	
									
	
								}
								elsif($size_diff <= 1)
								{
									warn "WARNING: Multiple Validation Results!\n";
									print "$newline\n";
									print "$this_pos\t$indel_type_454-$indel_size_454\t$ref\t$var\t$validation_status\n";																							
								}
							}							
						}

						

					}
				}
				
				
				## If no validation was made, try to refute indel ##
				
				if(!$newline)
				{
					## Check to see if indel was refuted ##
					my $validation_bases = my $non_indel_bases = 0;
					my $validation_lines = "";

					for(my $this_pos = $region_start - 1; $this_pos <= $region_stop + 1; $this_pos++)
					{
						$position_key = "$chromosome\t$this_pos";
						
						if($results454{$position_key})
						{
							my @results = split(/\t/, $results454{$position_key});
							
							my $ref = $results[2];
							my $var = $results[3];
							my $normal_gt = $results[7];
							my $tumor_gt = $results[11];
							my $validation_status = $results[12];						
	
							$validation_bases++;
							
							if($validation_status && $validation_status eq "Reference")
							{
								$non_indel_bases++;
								$validation_lines .= "\t$results454{$position_key}\n";
							}													
						}

					}

					if($validation_bases && $non_indel_bases == $validation_bases)
					{
#						print "$chromosome\t$position\t$var_type refuted by $non_indel_bases non-indel bases\n";
#						print $validation_lines . "\n";
						my $validation_status = "WildType";
						$lineContents[24] = "WildType";
						## Update the newline ##

						## Build a new line for MAF ##	
						$newline = "";
						my $numContents = @lineContents;
						for(my $colCounter = 0; $colCounter < $numContents; $colCounter++)
						{
							$lineContents[$colCounter] = "" if(!$lineContents[$colCounter]);
							$newline .= $lineContents[$colCounter] . "\t";
						}
						
						$stats{'changed from ' . $current_val_status . ' to ' . $validation_status}++;	
					}

				}
				
				
				if($newline)
				{
					print OUTFILE "$newline\n";
				}
				else
				{
					print OUTFILE "$line\n";
				}
				
			}
			else
			{
					## Variant not targeted - print original line ##
					print OUTFILE "$line\n";								
			}
	

		}

	}

	close($input);	
	
	print $stats{'num_targeted_positions'} . " positions targeted for validation\n";
	print $stats{'num_targeted_results'} . " had results in the status file\n";
	print $stats{'indels_with_results'} . " indels had results\n" if($stats{'indels_with_results'});
	print $stats{'num_type_matched'} . " matched the variant type\n" if($stats{'num_type_matched'});
	print $stats{'num_size_matched'} . " matched the variant size\n" if($stats{'num_size_matched'});

	## Print all of the changed stats ##

	foreach my $key (keys %stats)
	{
		if($key =~ "changed")
		{
			print $stats{$key} . " " . $key . "\n";
		}
	}
	
	close(OUTFILE);

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



#############################################################
# parse_file - takes input file and parses it
#
#############################################################

sub load_target_positions
{
	my $FileName = shift(@_);

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	my %results = ();
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		my @lineContents = split(/\t/, $line);			
		my $chromosome = $lineContents[0];
		my $position = $lineContents[1];
		my $position_key = "$chromosome\t$position";
		$results{$position_key} = 1;
	}

	close($input);
	
	return(%results);	
}




#############################################################
# parse_file - takes input file and parses it
#
#############################################################

sub load_454_results
{
	my $FileName = shift(@_);

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	my %results = ();
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		my @lineContents = split(/\t/, $line);			
		my $chromosome = $lineContents[0];
		my $position = $lineContents[1];
		my $ref_base = $lineContents[2];
		my $var_base = $lineContents[3];
		my $somatic_status = $lineContents[12];
		my $key = "$chromosome\t$position";#\t$ref_base\t$var_base";#\t$ref_base\t$var_base";
		my $result = $line;
		$results{$key} = $result;

		if($somatic_status eq "Reference")
		{
			my $position_key = "$chromosome\t$position";
			$results{$position_key} = $result;
		}
	}

	close($input);
	
	return(%results);	
}



1;

