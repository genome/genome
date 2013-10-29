
package Genome::Model::Tools::SnpArray::SegmentsToExons;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# SegmentsToExons - merges adjoining segments of similar copy number; distinguishes amplifications and deletions
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

my %stats = ();

class Genome::Model::Tools::SnpArray::SegmentsToExons {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		segments_file	=> { is => 'Text', doc => "Segments with p-values from running CBS on data", is_optional => 0, is_input => 1 },		
		exons_file	=> { is => 'Text', doc => "BED file of exons in chrom start stop gene format", is_optional => 0, is_input => 1 },		
		amp_threshold	=> { is => 'Text', doc => "Minimum seg_mean threshold for amplification", is_optional => 1, is_input => 1, default => 0.25},
		del_threshold	=> { is => 'Text', doc => "Maximum seg_mean threshold for deletion", is_optional => 1, is_input => 1, default => -0.25},
		size_threshold	=> { is => 'Text', doc => "Fraction of chromosome arm length above which an event is considered large-scale", is_input => 1, default => 0.25},
		output_file	=> { is => 'Text', doc => "Base name for output", is_optional => 1, is_input => 1},
		verbose	=> { is => 'Text', doc => "If set to 1, use for verbose output", is_optional => 1, is_input => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Merges adjoining segments of similar copy number"                 
}

sub help_synopsis {
    return <<EOS
This command assigns a copy number status to each exon in a BED file
EXAMPLE:	gmt snp-array segments-to-exons --segments-file snpArray.cbs.segments.tsv --output-basename snpArray.cbs.segments.exons
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

sub execute
{                               # replace with real execution logic.
	my $self = shift;

	## Get required parameters ##
	my $segments_file = $self->segments_file;
	my $exons_file = $self->exons_file;
	my $output_file = $self->output_file;

	## Get thresholds ##
	
	my $amp_threshold = $self->amp_threshold;
	my $del_threshold = $self->del_threshold;
	my $size_threshold = $self->size_threshold;

	## Reset various statistic counters ##
	my %stats = ();


	## Load exons ##
	warn "Loading exon targets...\n";
	my %exon_targets = load_exons($exons_file);

	## Parse the segments file ##
	print "Loading copy number segments...\n";
	my %segments = load_segments($segments_file);


	## Save exons with results ##

	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";	

	my %results = ();
	
	foreach my $chrom (sort byChrom keys %exon_targets)
	{
		warn "chr$chrom\n";
		my @exon_targets = split(/\n/, $exon_targets{$chrom});

		my %exons_with_results = ();
		
		foreach my $exon_target (@exon_targets)
		{
			$stats{'num_exons'}++;
			
			my ($target_chrom, $target_chr_start, $target_chr_stop, $target_name) = split(/\t/, $exon_target);
	
			if($segments{$chrom})# && $stats{'num_exons'} < 200)
			{
				my $result = find_match($target_chrom, $target_chr_start, $target_chr_stop, $target_name, $segments{$chrom}, $self);
				
				if($result)
				{
					$stats{'num_exons_with_results'}++;
#					print join("\t", $exon_target, $result) . "\n";
					
					$exons_with_results{$target_name} .= "\n" if($exons_with_results{$target_name});
					$exons_with_results{$target_name} .= join("\t", $target_chr_start, $target_chr_stop, $result);
				}
			}

		}

		## Go through each gene and grab its results ##
		
		foreach my $gene (keys %exons_with_results)
		{
			$stats{'num_genes_with_results'}++;
			my $gene_chr_start = my $gene_chr_stop = my $num_exons = my $seg_mean_sum = 0;
			my $gene_amp_bp = my $gene_del_bp = my $gene_neut_bp = 0;
			my $gene_amp_count = my $gene_del_count = my $gene_neut_count = 0;
			
			my @exonsWithResults = split(/\n/, $exons_with_results{$gene});
			
			foreach my $exon_result (@exonsWithResults)
			{
				my ($exon_chr_start, $exon_chr_stop, $event_type, $seg_mean) = split(/\t/, $exon_result);
				my $exon_size = $exon_chr_stop - $exon_chr_start + 1;
				
				## Save Start/STop ##
				$gene_chr_start = $exon_chr_start if(!$gene_chr_start || $exon_chr_start < $gene_chr_start);
				$gene_chr_stop = $exon_chr_stop if(!$gene_chr_stop || $exon_chr_stop > $gene_chr_stop);
				
				if($event_type eq "amp")
				{
					$gene_amp_bp += $exon_size;
					$gene_amp_count++;
				}
				elsif($event_type eq "del")
				{
					$gene_del_bp += $exon_size;
					$gene_del_count++;				
				}
				else
				{
					$gene_neut_bp += $exon_size;
					$gene_neut_count++;
				}

				## Save Seg_mean Sum ##
				$seg_mean_sum += ($seg_mean * $exon_size);

			}

			## Calculate average seg mean ##
			
			my $avg_seg_mean = $seg_mean_sum / ($gene_amp_bp + $gene_neut_bp + $gene_del_bp);
			$avg_seg_mean = sprintf("%.3f", $avg_seg_mean);

			## Determine if most exons were amp, del, or whatevs ##
			
			my @exon_results = ();
			my $resultCounter = 0;
			
			if($gene_amp_count)
			{
				$exon_results[$resultCounter] = join("\t", "amp", $gene_amp_count, $gene_amp_bp);
				$resultCounter++;
			}
			
			if($gene_del_count)
			{
				$exon_results[$resultCounter] = join("\t", "del", $gene_del_count, $gene_del_bp);				
				$resultCounter++;
			}
			
			if($gene_neut_count)
			{
				$exon_results[$resultCounter] = join("\t", "neut", $gene_neut_count, $gene_neut_bp);				
				$resultCounter++;				
			}


			@exon_results = sort byExonsThenBases @exon_results;
			
			sub byExonsThenBases
			{
				my ($type_a, $exons_a, $bases_a) = split(/\t/, $a);
				my ($type_b, $exons_b, $bases_b) = split(/\t/, $b);
				$exons_b <=> $exons_a
				or
				$bases_b <=> $bases_a;
			}

			## Get the predominant result ##
			
			my ($cns_type, $cns_exons, $cns_bases) = split(/\t/, $exon_results[0]);
			$results{$cns_type}++;

			print OUTFILE join("\t", $chrom, $gene_chr_start, $gene_chr_stop, $gene, $avg_seg_mean, $cns_type, $cns_exons, $cns_bases) . "\n";
		}


	}
	
	close(OUTFILE);
	
	print $stats{'num_exons'} . " exons in BED file\n";
	print $stats{'num_exons_with_results'} . " were matched to a copynumber result\n";
	print $stats{'num_genes_with_results'} . " genes with copynumber result\n";

	## Print summary of results ##
	
	foreach my $result (sort keys %results)
	{
		print $results{$result} . "\t" . $result . "\n";
	}
	
	
	sub byChrom
	{
		(my $chrom_a) = split(/\t/, $a);
		(my $chrom_b) = split(/\t/, $b);
		$chrom_a = 23 if($a eq "X");
		$chrom_a = 24 if($a eq "Y");
		$chrom_a = 25 if($a eq "MT");
		$chrom_a =~ s/[^0-9]//g;
	
		$chrom_b = 23 if($b eq "X");
		$chrom_b = 24 if($b eq "Y");
		$chrom_b = 25 if($b eq "MT");
		$chrom_b =~ s/[^0-9]//g;
	
		$chrom_a <=> $chrom_b;
	}



}



#############################################################
# parse_array - parses the file
#
#############################################################

sub find_match
{
	my ($chrom, $chr_start, $chr_stop, $type, $events, $self) = @_;
	
	my $size = $chr_stop - $chr_start + 1;
	my $best_overlap = my $best_event_type = my $best_event_mean = "";
	
	if($events)
	{
		my @events = split(/\n/, $events);		
		my $overlap = 0;

		foreach my $check_event (@events)
		{
			my ($event_chr_start, $event_chr_stop, $event_seg_mean) = split(/\t/, $check_event);
			my $event_size = $event_chr_stop - $event_chr_start + 1;

			## Determine event type ##
			
			my $event_type = "neut";
			
			if($event_seg_mean >= $self->amp_threshold)
			{
				$event_type = "amp";
			}
			elsif($event_seg_mean <= $self->del_threshold)
			{
				$event_type = "del";
			}

			if(1)
			{
				if($event_chr_stop >= $chr_start && $event_chr_start <= $chr_stop)
				{
					## Determine the overlap ##

					if($event_chr_start >= $chr_start && $event_chr_stop <= $chr_stop)
					{
						## Check event completely within sought event #
						$overlap = $event_chr_stop - $event_chr_start + 1;
					}
					elsif($chr_start >= $event_chr_start && $chr_stop <= $event_chr_stop)
					{
						## Sought event completely within check event #
						$overlap = $chr_stop - $chr_start + 1;
					}
					elsif($event_chr_start < $chr_start && $event_chr_stop >= $chr_start)
					{
						## Check event overlaps at front edge ##
						$overlap = $event_chr_stop - $chr_start;
					}
					elsif($event_chr_stop > $chr_stop && $event_chr_start <= $chr_stop)
					{
						## Check event overlaps at back edge ##
						$overlap = $chr_stop - $event_chr_start;
					}

					if(!$best_overlap || $overlap >= $best_overlap)
					{
						$best_event_type = $event_type;
						$best_overlap = $overlap;
						$best_event_mean = $event_seg_mean;
					}
				}
			}
		}
		

		
	}
	
#	return(join("\t", $best_event_type, $best_overlap));
	my $best_result = "";
	$best_result = "$best_event_type\t$best_event_mean" if($best_event_type);
	return($best_result);
	
}



################################################################################################
# Execute - the main program logic
#
################################################################################################

sub load_segments
{
	my %segments_by_chrom = ();
	
	my $segments_file = shift(@_);
	my $input = new FileHandle ($segments_file);
	my $lineCounter = 0;

	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($id, $sample, $chrom, $chr_start, $chr_stop, $num_mark, $seg_mean, $bstat, $p_value, $lcl, $ucl) = split(/\t/, $line);

		if($id eq "ID")
		{
			## Skip header ##
		}
		else
		{
			$chrom = "X" if($chrom eq "23");
			$chrom = "Y" if($chrom eq "24");
			$chrom =~ s/\"//g;		
			
			$segments_by_chrom{$chrom} .= "\n" if($segments_by_chrom{$chrom});
			$segments_by_chrom{$chrom} .= join("\t", $chr_start, $chr_stop, $seg_mean);
		}
	}
	
	close($input);
	
	return(%segments_by_chrom);

}




#############################################################
# parse_array - parses the file
#
#############################################################

sub load_exons
{
	my %events_by_chrom = ();
	my $FileName = shift(@_);

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	my $region_size_sum = my $region_size_num = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
		
		my ($chrom, $chr_start, $chr_stop, $gene) = split(/\t/, $line);

#		my $name = join(":", $chrom, $chr_start, $chr_stop, $gene);
		
		$events_by_chrom{$chrom} .= "\n" if($events_by_chrom{$chrom});
#		$events_by_chrom{$chrom} .= join("\t", $chrom, $chr_start, $chr_stop, $name);
		$events_by_chrom{$chrom} .= join("\t", $chrom, $chr_start, $chr_stop, $gene);

	}
	close($input);

	return(%events_by_chrom);
}


###############################################################################
# commify - add appropriate commas to long integers
###############################################################################

sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}







return(1);

