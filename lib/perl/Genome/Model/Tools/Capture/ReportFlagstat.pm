
package Genome::Model::Tools::Capture::ReportFlagstat;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ReportFlagstat - Compare tumor versus normal models to find somatic events
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	12/09/2009 by D.K.
#	MODIFIED:	12/09/2009 by D.K.
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

class Genome::Model::Tools::Capture::ReportFlagstat {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		file_list	=> { is => 'Text', doc => "List of all Flagstat files" , is_optional => 0},
		output_file	=> { is => 'Text', doc => "Output file for summary report" , is_optional => 1},
		read_length	=> { is => 'Text', doc => "Read length for total bp estimate", is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Summarizes BAM flag statistics across multiple files"                 
}

sub help_synopsis {
    return <<EOS
Summarizes BAM flag statistics across multiple files
EXAMPLE:	gmt capture report-flagstat ...
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
	my $file_list = $self->file_list;

	## Reset statistics ##
	$stats{'num_files'} = 0;

	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open output file: $!\n";
		
		if($self->read_length)
		{
			print OUTFILE "bam_file\tnum_reads\tpct_paired\tpct_frag\tnum_mapped\tpct_map\tnum_dup\tpct_dup\test_GBP\n";						
		}
		else
		{
			print OUTFILE "bam_file\tnum_reads\tpct_paired\tpct_frag\tnum_mapped\tpct_map\tnum_dup\tpct_dup\n";			
		}

	}

#	print "sample_name\ttotal_regions\tcov_20x\tcov_10x\tcov_8x\tcov_1x\tcov_0x\n";

	my $total_files = my $total_reads = my $total_mapped = my $total_dups = my $total_pair = my $total_frag = 0;

	my $input = new FileHandle ($file_list);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my @lineContents = split(/\t/, $line);
		my $numContents = @lineContents;
		
		my $flagstat_file = "";
		my $sample_name = "";
		
		if($numContents == 2)
		{
			($sample_name, $flagstat_file) = split(/\t/, $line);
		}
		else
		{
			($flagstat_file) = split(/\t/, $line);
		}

		$stats{'num_files'}++;
		
		my @pathContents = split(/\//, $flagstat_file);
		$numContents = @pathContents;
		my $local_filename = $pathContents[$numContents - 1];
		$local_filename = $sample_name if($sample_name);
		
		if(-e $flagstat_file)
		{
			my %bam_stats = parse_flagstat_file($flagstat_file);
			
			my $num_reads = $bam_stats{'in total'};
			my $qc_failure = $bam_stats{'QC failure'};
			my $num_dups = $bam_stats{'duplicates'};
			my $num_mapped = $bam_stats{'mapped'};
			my $num_paired = $bam_stats{'paired in sequencing'};
			my $pair_mapped = $bam_stats{'with itself and mate mapped'};
			my $pair_proper = $bam_stats{'properly paired'};
			my $pair_singleton = $bam_stats{'singletons'};
	
			$num_reads = 0 if(!$num_reads);
			$qc_failure = 0 if(!$qc_failure);
			$num_dups = 0 if(!$num_dups);
			$num_mapped = 0 if(!$num_mapped);
			$num_paired = 0 if(!$num_paired);
			
			$pair_mapped = 0 if(!$pair_mapped);
			$pair_proper = 0 if(!$pair_proper);
			$pair_singleton = 0 if(!$pair_singleton);
			my $num_fragment = $num_reads - $num_paired;
			
			my $map_rate = my $dup_rate = "0.00%";
			my $pct_fragment = my $pct_paired = my $pct_proper = "0.00%";
			
			if($num_reads)
			{
				$map_rate = sprintf("%.2f", ($num_mapped / $num_reads * 100)) . '%';
				$dup_rate = sprintf("%.2f", ($num_dups / $num_reads * 100)) . '%';
				$pct_fragment = sprintf("%.2f", ($num_fragment / $num_reads * 100)) . '%';
				$pct_paired = sprintf("%.2f", ($num_paired / $num_reads * 100)) . '%';
				$pct_proper = sprintf("%.2f", ($pair_proper / $num_paired * 100)) . '%' if($num_paired > 0);
			
				$total_reads += $num_reads;
				$total_mapped += $num_mapped;
				$total_dups += $num_dups;
				$total_pair += $num_paired;
				$total_frag += $num_fragment;
				$total_files++;
			}
	
	
	
			$num_dups = commify($num_dups);
			$num_paired = commify($num_paired);
			$num_fragment = commify($num_fragment);
			$num_mapped = commify($num_mapped);
	
			if($self->read_length)
			{
				my $est_gigabases = sprintf("%.2f", ($num_reads * $self->read_length / 1000000000));
				$num_reads = commify($num_reads);
				print "$local_filename\t$num_reads\t$pct_paired\t$pct_fragment\t$num_mapped\t$map_rate\t$num_dups\t$dup_rate\t$pct_proper\t$est_gigabases\n";						
				print OUTFILE "$local_filename\t$num_reads\t$pct_paired\t$pct_fragment\t$num_mapped\t$map_rate\t$num_dups\t$dup_rate\t$est_gigabases\n" if($self->output_file);
			}
			else
			{
				$num_reads = commify($num_reads);
				print "$local_filename\t$num_reads\t$pct_paired\t$pct_fragment\t$num_mapped\t$map_rate\t$num_dups\t$dup_rate\t$pct_proper\n";
				print OUTFILE "$local_filename\t$num_reads\t$pct_paired\t$pct_fragment\t$num_mapped\t$map_rate\t$num_dups\t$dup_rate\n" if($self->output_file);
			}

		}
	}

	close($input);

	
	print $stats{'num_files'} . " flagstat files\n";

	## Print some summary statistics ##
	
	my $pct_map = sprintf("%.2f", ($total_mapped / $total_reads * 100)) . '%';
	my $pct_dup = sprintf("%.2f", ($total_dups / $total_reads * 100)) . '%';
	my $pct_pair = sprintf("%.2f", ($total_pair / $total_reads * 100)) . '%';
	my $pct_frag = sprintf("%.2f", ($total_frag / $total_reads * 100)) . '%';	


	my $avg_reads = sprintf("%d", $total_reads / $total_files);
	my $est_gigabases;
	
	if($self->read_length)
	{
		$est_gigabases = sprintf("%.2f", ($avg_reads * $self->read_length / 1000000000));
	}

	$avg_reads = commify($avg_reads);

	print "Avg. Reads Per File:\t" . $avg_reads . "\n";
	print "Pct. Paired-end:\t" . $pct_pair . "\n";
	print "Pct. Fragment-end:\t" . $pct_frag . "\n";
	print "Avg. Percent Mapped:\t" . $pct_map . "\n";
	print "Avg. Duplication:\t" . $pct_dup . "\n";
	if($self->read_length)
	{
		print "Avg. Gigabases Per:\t" . $est_gigabases . "\n";
	}

	if($self->output_file)
	{
		print OUTFILE "Avg. Reads Per File:\t" . $avg_reads . "\n";
		print OUTFILE "Pct. Paired-end:\t" . $pct_pair . "\n";
		print OUTFILE "Pct. Fragment-end:\t" . $pct_frag . "\n";
		print OUTFILE "Avg. Percent Mapped:\t" . $pct_map . "\n";
		print OUTFILE "Avg. Duplication:\t" . $pct_dup . "\n";
		if($self->read_length)
		{
			print OUTFILE "Avg. Gigabases:\t" . $est_gigabases . "\n";
		}
	}
	
	close(OUTFILE) if($self->output_file);


	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




#############################################################
# ParseFile - takes input file and parses it
#
#############################################################

sub parse_flagstat_file
{
	(my $stats_file) = @_;
	
	## Parse the file ##

	my $input = new FileHandle ($stats_file);
	my $lineCounter = 0;
	
	my %cov_stats = ();
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		(my $num_reads) = split(/\s+/, $line);

		my $category = $line;
		$category =~ s/$num_reads\s//;

		if($category =~ '\+')
		{
			$category = "in total" if($category =~ 'in total');
			$category = "duplicates" if($category =~ 'duplicates');
			$category = "mapped" if($category =~ 'mapped' && !($category =~ 'mapped to'));
			$category = "paired in sequencing" if($category =~ 'paired in');
			$category = "read1" if($category =~ 'read1');
			$category = "read2" if($category =~ 'read2');
			$category = "properly paired" if($category =~ 'properly');			
		}
		
		## Remove stuff with parentheses ##
		my $split_char = " \\(";
		($category) = split(/$split_char/, $category);
		
		$cov_stats{$category} = $num_reads if($category);
	}
	
	close($input);

	return(%cov_stats);
	
}


sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}

1;

