
package Genome::Model::Tools::SnpArray::CompareCopySegments;     # rename this when you give the module file a different name <--

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

class Genome::Model::Tools::SnpArray::CompareCopySegments {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		array_events	=> { is => 'Text', doc => "List of merged CNA events from SNP array data, e.g. snpArray.cbs.segments.events.tsv", is_optional => 0, is_input => 1 },
		sequence_events	=> { is => 'Text', doc => "List of merged CNA events from SNP array data, e.g. varScan.output.copynumber.cbs.segments.events.tsv", is_optional => 0, is_input => 1 },
		wgs_events	=> { is => 'Text', doc => "List of merged CNA events from WGS data, e.g. cnvHmm.cbs.segments.merged.events.tsv", is_optional => 1, is_input => 1 },
		event_size	=> { is => 'Text', doc => "Either large-scale or focal or all", is_optional => 1, is_input => 1, default => "all"},
		output_file	=> { is => 'Text', doc => "Output file for comparison result", is_optional => 1, is_input => 1},
		output_hits	=> { is => 'Text', doc => "Output file shared events", is_optional => 1, is_input => 1},
		output_misses	=> { is => 'Text', doc => "Output file for events from only one source", is_optional => 1, is_input => 1},
		verbose	=> { is => 'Text', doc => "Prints verbosely if set to 1", is_optional => 1, is_input => 1}
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Compares copy number alterations detected by SNP array versus sequence data"                 
}

sub help_synopsis {
    return <<EOS
This command compares copy number alterations detected by SNP array versus sequence data
EXAMPLE:	gmt snp-array compare-copy-segments --array-events snpArray.cbs.segments.events.tsv --sequence-events varScan.output.copynumber.cbs.segments.events.tsv
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
	my $array_events_file = $self->array_events;
	my $sequence_events_file = $self->sequence_events;
	my $output_file = $self->output_file;

	my %stats = ();
	$stats{'num_sequence_events'} = $stats{'num_supported_by_array'} = $stats{'num_sequence_amps'} = $stats{'num_amps_supported_by_array'} = $stats{'num_sequence_dels'} = $stats{'num_dels_supported_by_array'} = $stats{'num_array_events'} = $stats{'num_supported_by_sequence'} = $stats{'num_array_amps'} = $stats{'num_amps_supported_by_sequence'} = $stats{'num_array_dels'} = $stats{'num_dels_supported_by_sequence'} = 0;		

	print "Loading SNP array events...\n" if($self->verbose);
	my %array_events = load_events($array_events_file, $self);

	print "Loading sequence events...\n" if($self->verbose);
	my %sequence_events = load_events($sequence_events_file, $self);

	my %wgs_events = load_events($self->wgs_events, $self) if($self->wgs_events);

	## Open hits and misses files if specified ##
	
	if($self->output_hits)
	{
		open(HITS, ">" . $self->output_hits) or die "Can't open outfile: $!\n";
	}

	if($self->output_misses)
	{
		open(MISSES, ">" . $self->output_misses) or die "Can't open outfile: $!\n";
	}

	if($self->output_file)
	{
		open(OUTFILE, ">" . $output_file) or die "Can't open outfile: $!\n";
	}


	## Go through the events by chromosome ##
	if($self->wgs_events)
	{
		## Do A three way comparison ##
		foreach my $chrom (sort keys %array_events)
		{
			if($array_events{$chrom} && $sequence_events{$chrom} && $wgs_events{$chrom})
			{
				my %chrom_stats = three_way_comparison($sequence_events{$chrom}, $array_events{$chrom}, $wgs_events{$chrom}, 0);

				foreach my $key (keys %chrom_stats)
				{
					if($key ne 'output')
					{
						$stats{$key} += $chrom_stats{$key};						
					}
					else
					{
#						print "Printing to outfile: " . $chrom_stats{$key} . "\n";
						print OUTFILE $chrom_stats{$key} . "\n";
					}

				}
			}
		}
		
		print "events\tsupp_both\tarray_only\twgs_only\tnot_supp\n";
		$stats{'num_events'} = 0 if(!$stats{'num_events'});
		$stats{'supported_both'} = 0 if(!$stats{'supported_both'});
		$stats{'supported_array'} = 0 if(!$stats{'supported_array'});
		$stats{'supported_wgs'} = 0 if(!$stats{'supported_wgs'});
		$stats{'not_supported'} = 0 if(!$stats{'not_supported'});
		print join("\t", $stats{'num_events'}, $stats{'supported_both'}, $stats{'supported_array'}, $stats{'supported_wgs'}, $stats{'not_supported'}) . "\n";

	}
	else
	{
		## Do a two-way comparison ##
		foreach my $chrom (sort keys %array_events)
		{
			if($array_events{$chrom} && $sequence_events{$chrom})
			{
				## Comparison 1: Specificity: How many Sequence-based Events are Supported by Array Calls? ##
				## The third parameter (1) tells it to print overlaps ##
				my %chrom_stats = compare_events($array_events{$chrom}, $sequence_events{$chrom}, 0);
	
				## Update the Overall Totals ##
				$stats{'num_sequence_events'} += $chrom_stats{'num_events'};
				$stats{'num_sequence_amps'} += $chrom_stats{'num_amps'};
				$stats{'num_sequence_dels'} += $chrom_stats{'num_dels'};
	
				$stats{'num_supported_by_array'} += $chrom_stats{'num_supported'};
				$stats{'num_amps_supported_by_array'} += $chrom_stats{'num_amps_supported'};
				$stats{'num_dels_supported_by_array'} += $chrom_stats{'num_dels_supported'};
	
	#			print join("\t", $chrom, $chrom_stats{'num_events'}, $chrom_stats{'num_supported'}) . "\n";
	
				if($self->output_file)
				{
					print OUTFILE $stats{'output'} . "\n";					
				}

	
				## Print hits and misses ##
			
				if($self->output_hits)
				{
					print HITS "$chrom_stats{'hits'}";
				}
			
				if($self->output_misses)
				{
					print MISSES "$chrom_stats{'misses'}";
				}
				
	
				## Comparison 2: Sensitivity: How many Array-based Events are Detected by Sequence calls? ##
				%chrom_stats = ();
				%chrom_stats = compare_events($sequence_events{$chrom}, $array_events{$chrom}, 0);
				$stats{'num_array_events'} += $chrom_stats{'num_events'};
				$stats{'num_array_amps'} += $chrom_stats{'num_amps'};
				$stats{'num_array_dels'} += $chrom_stats{'num_dels'};			
	
				$stats{'num_supported_by_sequence'} += $chrom_stats{'num_supported'};
				$stats{'num_amps_supported_by_sequence'} += $chrom_stats{'num_amps_supported'};
				$stats{'num_dels_supported_by_sequence'} += $chrom_stats{'num_dels_supported'};			
	
	
	
			}
		}		
	
		print $stats{'num_supported_by_array'} . " of " . $stats{'num_sequence_events'} . " sequence events supported by array\n";	
		print $stats{'num_amps_supported_by_array'} . " of " . $stats{'num_sequence_amps'} . " sequence amplifications supported by array\n";	
		print $stats{'num_dels_supported_by_array'} . " of " . $stats{'num_sequence_dels'} . " sequence deletions supported by array\n";	
	
		print "\n";
	
		print $stats{'num_supported_by_sequence'} . " of " . $stats{'num_array_events'} . " array events supported by sequence\n";	
		print $stats{'num_amps_supported_by_sequence'} . " of " . $stats{'num_array_amps'} . " array amplifications supported by sequence\n";	
		print $stats{'num_dels_supported_by_sequence'} . " of " . $stats{'num_array_dels'} . " array deletions supported by sequence\n";	

	
		if($self->output_file)
		{
			print OUTFILE "sequence_events\tsupported_by_array\tamps\tsupported\tdels\tsupported\tarray_events\tdetected_by_sequence\tamps\tdetected\tdels\tdetected\n";
	#		print OUTFILE join("\t", $stats{'num_sequence_events'}, $stats{'num_supported_by_array'}, $stats{'num_sequence_amps'}, $stats{'num_amps_supported_by_array'}, $stats{'num_sequence_dels'}, $stats{'num_dels_supported_by_array'});
	#		print OUTFILE "\t";
	#		print OUTFILE join("\t", $stats{'num_array_events'}, $stats{'num_supported_by_sequence'}, $stats{'num_array_amps'}, $stats{'num_amps_supported_by_sequence'}, $stats{'num_array_dels'}, $stats{'num_dels_supported_by_sequence'});
	#		print OUTFILE "\n";
		}


	}



	if($self->output_hits)
	{
		close(HITS);
	}

	if($self->output_misses)
	{
		close(MISSES);
	}

	if($self->output_file)
	{
		close(OUTFILE);
	}

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


################################################################################################
# Load Genotypes
#
################################################################################################

sub load_events
{                               # replace with real execution logic.
	my $events_file = shift(@_);
	my $self = shift(@_);
	my %events = ();

	my %type_counts = ();
	
	my $input = new FileHandle ($events_file);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		if($lineCounter > 1)
		{
#			my ($id, $sample, $chrom, $chr_start, $chr_stop, $num_mark, $seg_mean, $bstat, $p_value, $lcl, $ucl, $event_size, $event_type) = split(/\t/, $line);
#chrom\tchr_start\tchr_stop\tseg_mean\tnum_segments\tnum_markers\tp_value\tevent_type\tevent_size\tsize_class\tchrom_arm\tarm_fraction\tchrom_fraction
#			my ($chrom, $chr_start, $chr_stop, $seg_mean, $num_segments, $num_markers, $p_value, $event_type, $event_size_bp, $event_size, $chrom_arm) = split(/\t/, $line);
			
			my ($chrom, $chr_start, $chr_stop, $seg_mean, $num_segments, $num_mark, $p_value, $event_type, $event_bases, $event_size, $chrom_arm) = split(/\t/, $line);

			if($self->event_size eq 'all' || $event_size eq $self->event_size)
			{
				$events{$chrom} .= "\n" if($events{$chrom});
				$events{$chrom} .= join("\t", $chrom, $chr_start, $chr_stop, $num_mark, $seg_mean, $p_value, $event_size, $event_type);
	
				$type_counts{"$event_size $event_type"}++;				
			}

		}

	}
	close($input);


	foreach my $event_type (sort keys %type_counts)
	{
		print "$type_counts{$event_type} $event_type, " if($self->verbose);
	}
	
	print "\n" if($self->verbose);
	
	return(%events);                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




################################################################################################
# Load Genotypes
#
################################################################################################

sub compare_events
{
	my ($array_events, $sequence_events, $print_overlaps) = @_;
	
	my %chrom_stats = ();
	$chrom_stats{'num_events'} = $chrom_stats{'num_supported'} = 0;
	$chrom_stats{'num_amps'} = $chrom_stats{'num_amps_supported'} = 0;
	$chrom_stats{'num_dels'} = $chrom_stats{'num_dels_supported'} = 0;
	
	my @array_events = split(/\n/, $array_events);
	my @sequence_events = split(/\n/, $sequence_events);
	
	foreach my $sequence_event (@sequence_events)
	{
		my ($chrom, $chr_start, $chr_stop, $num_mark, $seg_mean, $p_value, $event_size, $event_type) = split(/\t/, $sequence_event);
		
		$chrom_stats{'num_events'}++;

		if($event_type =~ 'amp')
		{
			$chrom_stats{'num_amps'}++;
		}
		elsif($event_type =~ 'del')
		{
			$chrom_stats{'num_dels'}++;
		}
		
		## Look for supporting events on the array ##
		
		my $array_supported_flag = 0;		
		my $array_overlaps = "";
		
		my $best_overlap_bp = my $best_overlap = 0;
		
		foreach my $array_event (@array_events)
		{
			my ($array_chrom, $array_chr_start, $array_chr_stop, $array_num_mark, $array_seg_mean, $array_p_value, $array_event_size, $array_event_type) = split(/\t/, $array_event);
			
			
			
			## Match Chromosome ##
			if($array_chrom eq $chrom)
			{
				## Check Positional Overlap ##
				if($array_chr_stop >= $chr_start && $array_chr_start <= $chr_stop)
				{
					my $overlap = get_overlap($chr_start, $chr_stop, $array_chr_start, $array_chr_stop);

					if($overlap > $best_overlap_bp)
					{
						$best_overlap = join("\t", "", $array_chrom, $array_chr_start, $array_chr_stop, $array_seg_mean, $array_event_size, $array_event_type);
						$best_overlap_bp = $overlap;
					}
				}
			}
		}
		
		if($best_overlap)
		{
			$chrom_stats{'num_supported'}++;

			if($event_type =~ 'amp')
			{
				$chrom_stats{'num_amps_supported'}++;
			}
			elsif($event_type =~ 'del')
			{
				$chrom_stats{'num_dels_supported'}++;
			}

#			if($print_overlaps)
#			{
				print join("\t", $chrom, $chr_start, $chr_stop, $seg_mean, $event_type, $best_overlap, $best_overlap_bp) . "\n";				
#			}

			## Save hits ##
			$chrom_stats{'hits'} .= "\n" if($chrom_stats{'hits'});
			$chrom_stats{'hits'} .= join("\t", "EVENT", $sequence_event) . "\n";
			$chrom_stats{'hits'} .= $array_overlaps;
		}
		else
		{
			## Save misses ##
			$chrom_stats{'misses'} .= "\n" if($chrom_stats{'misses'});
			$chrom_stats{'misses'} .= join("\t", "EVENT", $sequence_event) . "\n";
		}
	}
	
	return(%chrom_stats);
}




################################################################################################
# Load Genotypes
#
################################################################################################

sub three_way_comparison
{
	my ($sequence_events, $array_events, $wgs_events, $print_overlaps) = @_;
	
	my %chrom_stats = ();
	$chrom_stats{'num_events'} = $chrom_stats{'num_supported'} = 0;
	$chrom_stats{'num_amps'} = $chrom_stats{'num_amps_supported'} = 0;
	$chrom_stats{'num_dels'} = $chrom_stats{'num_dels_supported'} = 0;
	
	my @sequence_events = split(/\n/, $sequence_events);
	my @array_events = split(/\n/, $array_events);
	my @wgs_events = split(/\n/, $wgs_events);
	
	foreach my $sequence_event (@sequence_events)
	{
		my ($chrom, $chr_start, $chr_stop, $num_mark, $seg_mean, $p_value, $event_size, $event_type) = split(/\t/, $sequence_event);
		
		$chrom_stats{'num_events'}++;

		if($event_type =~ 'amp')
		{
			$chrom_stats{'num_amps'}++;
		}
		elsif($event_type =~ 'del')
		{
			$chrom_stats{'num_dels'}++;
		}
		
		## Check for support on array ##
		
		my $best_array_overlap = get_best_overlap($chrom, $chr_start, $chr_stop, $event_size, $event_type, $array_events);
		my $best_wgs_overlap = get_best_overlap($chrom, $chr_start, $chr_stop, $event_size, $event_type, $wgs_events);
		my $comparison_result = "unknown";
		my $comparison_details = "";

		if($best_array_overlap && $best_wgs_overlap)
		{
			$comparison_result = "supported_both";
			$comparison_details = $best_array_overlap . "\t" . $best_wgs_overlap;
		}
		elsif($best_array_overlap)
		{
			$comparison_result = "supported_array";
			$comparison_details = $best_array_overlap;
		}
		elsif($best_wgs_overlap)
		{
			$comparison_result = "supported_wgs";
			$comparison_details = $best_wgs_overlap;
		}
		else
		{
			$comparison_result = "not_supported";
		}
		
		$chrom_stats{$comparison_result}++;
		
		$chrom_stats{'output'} .= "\n" if($chrom_stats{'output'});
		$chrom_stats{'output'} .= $sequence_event . "\t" . $comparison_result . "\t" . $comparison_details;
	}
	
	return(%chrom_stats);
}



################################################################################################
# Load Genotypes
#
################################################################################################

sub get_best_overlap
{
	my ($chrom, $chr_start, $chr_stop, $event_size, $event_type, $test_events) = @_;
	
	my @test_events = split(/\n/, $test_events);
	my $best_overlap_bp = my $best_overlap = 0;
	
	foreach my $test_event (@test_events)
	{
		my ($test_chrom, $test_chr_start, $test_chr_stop, $test_num_mark, $test_seg_mean, $test_p_value, $test_event_size, $test_event_type) = split(/\t/, $test_event);
		
		## Match Chromosome ##
		if($test_chrom eq $chrom && $test_event_type eq $event_type)
		{
			## Check Positional Overlap ##
			if($test_chr_stop >= $chr_start && $test_chr_start <= $chr_stop)
			{
				my $overlap = get_overlap($chr_start, $chr_stop, $test_chr_start, $test_chr_stop);

				if($overlap > $best_overlap_bp)
				{
					$best_overlap = join("\t", "", $test_chrom, $test_chr_start, $test_chr_stop, $test_seg_mean, $test_event_size, $test_event_type, $overlap);
					$best_overlap_bp = $overlap;
				}
			}
		}
	}
	
	return($best_overlap);
}


sub get_overlap
{
	my $overlap = 0;
	my ($start1, $stop1, $start2, $stop2) = @_;
	
	## 1111111111111
	##    22222222
	if($start2 >= $start1 && $stop2 <= $stop1)
	{
		$overlap = $stop2 - $start2 + 1;	
	}
	##     111111
	##    22222222		
	elsif($start1 >= $start2 && $stop1 <= $stop2)
	{
		$overlap = $stop1 - $start1 + 1;
	}
	## 111111
	##    22222222
	elsif($start1 <= $start2 && $stop1 >= $start2)
	{
		$overlap = $stop1 - $start2 + 1;
	}
	##       1111111111111
	##    22222222
	elsif($start2 <= $start1 && $stop2 >= $start1)
	{
		$overlap = $stop2 - $start1 + 1;
	}
	##      1111111
	##    22222222		
	elsif($start1 <= $stop2 && $stop1 >= $stop2)
	{
		$overlap = $stop2 - $start1 + 1;
	}
	## 111111111
	##    22222222		
	elsif($start2 <= $stop1 && $stop2 >= $stop1)
	{
		$overlap = $stop1 - $start2 + 1;	
	}
	
	return($overlap);
}

sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}


1;

