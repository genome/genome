
package Genome::Model::Tools::Varscan::ProcessValidationIndels;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ProcessValidation - Report the results of validation 
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	10/21/2010 by D.K.
#	MODIFIED:	10/21/2010 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

my %target_indels = ();

class Genome::Model::Tools::Varscan::ProcessValidationIndels {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		validation_snp_file		=> { is => 'Text', doc => "Varscan output file for validation data", is_optional => 0 },
		validation_indel_file		=> { is => 'Text', doc => "Varscan calls passing strand-filter in validation BAM (recommended)", is_optional => 0 },
		variants_file 	=> { is => 'Text', doc => "File of variants to report on", is_optional => 0 },
		output_file 	=> { is => 'Text', doc => "Output file for validation results", is_optional => 0 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Processes and reports on validation status of a list of indels"                 
}

sub help_synopsis {
    return <<EOS
Processes and reports on validation status of a list of indels
EXAMPLE:	gmt capture process-validation ...
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
	my $validation_snp_file = $self->validation_snp_file;
	my $validation_indel_file = $self->validation_indel_file;
	my $variants_file = $self->variants_file;
	my $output_file = $self->output_file;
	

	## Load the target indels ##

	print "Loading target indels...\n";	
	%target_indels = load_target_indels($variants_file);


	## Load SNP and indel validation calls at/around this position ##

	warn "Loading validation indels...\n";
	my %validation_indels = load_validation_results($validation_indel_file);

	warn "Loading validation SNPs...\n";
	my %validation_snps = load_validation_results($validation_snp_file);


	## Parse the list of target variants, looking for ones with validation results ##


	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	print OUTFILE "chrom\tchr_start\tchr_stop\tref\tvar\tcalled\tfilter\tstatus\tv_ref\tv_var\tnorm_reads1\tnorm_reads2\tnorm_freq\tnorm_call\ttum_reads1\ttum_reads2\ttum_freq\ttum_call\tsomatic_status\tgermline_p\tsomatic_p\n";

	## Reset statistics ##
	
	my %stats = ();



	## Parse the variant file ##

	my $input = new FileHandle ($variants_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($chrom, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $line);
		
		$stats{'num_variants'}++;

		## Determine indel type and size ##
		my $indel_type = my $indel_size = "";
		if($ref eq "0" || $ref eq "-" || length($var) > 1)
		{
			$indel_type = "INS";
			$indel_size = length($var);
		}
		else
		{
			$indel_type = "DEL";
			$indel_size = length($ref);
		}
		
		## Look for indels in validation output ##
		
		my $indel_results = my $snp_results = "";
		
		for(my $position = ($chr_start - $indel_size - 2); $position <= ($chr_stop + $indel_size + 2); $position++)
		{
			my $key = join("\t", $chrom, $position);

			if($validation_indels{$key})
			{
				my @validationIndel = split(/\t/, $validation_indels{$key});
				my $val_indel_type = my $val_indel_size = "";
				
				if($validationIndel[3] && ($validationIndel[3] =~ 'INS' || $validationIndel[3] =~ 'DEL'))
				{
					($val_indel_type, $val_indel_size) = split(/\-/, $validationIndel[3]);	
				}
				else
				{
					my $plus_or_minus = substr($validationIndel[3], 0, 1);
					my $indel_bases = substr($validationIndel[3], 1, 99);
					$val_indel_size = length($indel_bases);
					if($plus_or_minus eq '+')
					{
						$val_indel_type = "INS";
					}
					else
					{
						$val_indel_type = "DEL";
					}
				}
				
				if($val_indel_type eq $indel_type && $val_indel_size == $indel_size)
				{
					## Get indel distance ##
					
					my $indel_distance = min_indel_distance($chr_start, $chr_stop, $position);
					
					$indel_results .= "\n" if($indel_results);
					$indel_results .= "$indel_distance\t" . $validation_indels{$key};
				}


			}

			if($validation_snps{$key})
			{
				$snp_results .= "\n" if($snp_results);
				$snp_results .= $validation_snps{$key};
			}
		}

		## Print the indel along with results ##
		
		my $validation_result = my $validation_status = "";
		
		if($indel_results)
		{
			my @indelResults = split(/\n/, $indel_results);

			## Get the top matching indel result ##
			@indelResults = sort byDistanceThenReadcount(@indelResults);
			my @bestContents = split(/\t/, $indelResults[0]);
			my $numContents = @bestContents;

			## Append columns 2-x to the validation result ##
			
			for(my $colCounter = 4; $colCounter < $numContents; $colCounter++)
			{
				$validation_result .= "\t" if($validation_result);
				$validation_result .= $bestContents[$colCounter];
			}
			
			## Count this one ##
			$stats{'with_validation_result'}++;			
			$validation_status = $bestContents[13];
			
			$stats{'validated_' . $validation_status}++;


		}
		elsif($snp_results)
		{
			## Get the average reference-supporting read count ##
			my $num_positions = 0;
			my $normal_reads1_sum = my $tumor_reads1_sum = 0;
			
			my @snpResults = split(/\n/, $snp_results);
			foreach my $result (@snpResults)
			{
				my @resultContents = split(/\t/, $result);
				my $normal_reads1 = $resultContents[4];
				my $tumor_reads1 = $resultContents[8];
				$normal_reads1_sum += $normal_reads1;
				$tumor_reads1_sum += $tumor_reads1;
				$num_positions++;
			}
			
			if($num_positions)
			{
				my $avg_depth_normal = sprintf("%d", $normal_reads1_sum / $num_positions);
				my $avg_depth_tumor = sprintf("%d", $tumor_reads1_sum / $num_positions);
				
				$validation_result = join("\t", "N", $avg_depth_normal, "0", "0%", $ref, $avg_depth_tumor, "0", "0%", $ref, "Refuted", "1", "1");
				$stats{'indel_not_present'}++;
				$validation_status = "Refuted";
			}
		}
		
		if($validation_result)
		{
			print OUTFILE join("\t", $chrom, $chr_start, $chr_stop, $ref, $var, $validation_status, $validation_result) . "\n";
		}
		else
		{
			$stats{'no_coverage'}++;
			print OUTFILE join("\t", $chrom, $chr_start, $chr_stop, $ref, $var, "NoCoverage") . "\n";
		}

	}
	
	close($input);
	

	print $stats{'num_variants'} . " variants in $variants_file\n";
	print $stats{'no_coverage'} . " had no coverage\n";
	print $stats{'indel_not_present'} . " showed no indel in normal or tumor and were refuted\n";
	print $stats{'with_validation_result'} . " had a matching variant in normal and/or tumor\n";

	foreach my $key (sort keys %stats)
	{
		print "\t" . $stats{$key} . " $key\n" if($key =~ 'validated_');
	}

#	print $stats{'with_filtered_results'} . " had post-filter validation results\n";
#	print $stats{'with_unfiltered_results'} . " had unfiltered validation results\n";


	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



################################################################################################
# Process results - filter variants by type and into high/low confidence
#
################################################################################################

sub byDistanceThenReadcount
{
	my @temp = split(/\t/, $a);
	my $dist_a = $temp[0];
	my $count_a = $temp[6];
	
	@temp = ();
	@temp = split(/\t/, $b);
	my $dist_b = $temp[0];
	my $count_b = $temp[6];

	$dist_a <=> $dist_b
	or
	$count_b <=> $count_a;
}

################################################################################################
# Process results - filter variants by type and into high/low confidence
#
################################################################################################

sub min_indel_distance
{
	my ($chr_start, $chr_stop, $position) = @_;
	
	my $dist1 = abs($chr_start - $position);
	my $dist2 = abs($chr_start - $position);
	
	return($dist2) if($dist2 < $dist1);
	return($dist1);
}


################################################################################################
# Process results - filter variants by type and into high/low confidence
#
################################################################################################

sub load_validation_results
{
	my $filename = shift(@_);
	
	my %results = ();
	
	my $input = new FileHandle ($filename);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		(my $chrom, my $position) = split(/\t/, $line);
		
		my $key = join("\t", $chrom, $position);
		
		if($target_indels{$key})
		{
			$results{$key} = $line;			
		}

	}
	
	close($input);	

	return(%results);
}




################################################################################################
# Process results - filter variants by type and into high/low confidence
#
################################################################################################

sub load_target_indels
{
	my $filename = shift(@_);
	
	my %results = ();
	
	my $input = new FileHandle ($filename);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($chrom, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $line);
		
		
		## Determine indel type and size ##
		my $indel_type = my $indel_size = "";
		if($ref eq "0" || $ref eq "-" || length($var) > 1)
		{
			$indel_type = "INS";
			$indel_size = length($var);
		}
		else
		{
			$indel_type = "DEL";
			$indel_size = length($ref);
		}

		## Save this indel at nearby 

		for(my $position = ($chr_start - $indel_size - 2); $position <= ($chr_stop + $indel_size + 2); $position++)
		{
			my $key = join("\t", $chrom, $position);
			$results{$key} = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
		}

	}
	
	close($input);	

	return(%results);
}

1;

