
package Genome::Model::Tools::Capture::ReportRefcov;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ReportRefcov - Compare tumor versus normal models to find somatic events
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	12/09/2009 by D.K.
#	MODIFIED:	12/09/2009 by D.K.
#
#	NOTES:	This tool parses RefCov summary files to generate project-level coverage reports.  
#
#
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

## Declare global statistics hash ##

my %stats = ();

class Genome::Model::Tools::Capture::ReportRefcov {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		sample_stats_files	=> { is => 'Text', doc => "Tab-delimited list of sample names and refcov summary files" , is_optional => 0},
		sample_alignsum_files	=> { is => 'Text', doc => "Tab-delimited list of sample names and alignment summary files" , is_optional => 1},
		output_file	=> { is => 'Text', doc => "Output file for summary report" , is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Summarizes RefCov coverage across capture samples"                 
}

sub help_synopsis {
    return <<EOS
Summarizes RefCov coverage across capture samples
EXAMPLE:	gmt capture report-coverage ...
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
	my $sample_stats_files = $self->sample_stats_files;
	my %alignment_summaries = ();
	
	if($self->sample_alignsum_files)
	{
		%alignment_summaries = load_alignment_summaries($self->sample_alignsum_files);
	}

	my $alignsum_header = "unique_on_target\tdup_on_target\tunique_off_target\tdup_off_target\ton_target_rate";

	## Reset statistics ##
	$stats{'num_samples'} = 0;

	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open output file: $!\n";
#		print OUTFILE "sample_name\ttotal_regions\tcov_20x\tcov_10x\tcov_8x\tcov_1x\tcov_0x\n";
	}

	my $input = new FileHandle ($sample_stats_files);
	my $lineCounter = 0;
	my $header_printed = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		(my $sample_name, my $refcov_file) = split(/\t/, $line);
		$stats{'num_samples'}++;

		if(-e $refcov_file)
		{
			$stats{'with_stats'}++;
			(my $stats_header, my $sample_stats) = split(/\n/, parse_stats_file($refcov_file));

			## Print header if not done ##
			
			if(!$header_printed)
			{
				print "Sample\t$stats_header";
				print "\t$alignsum_header" if($self->sample_alignsum_files);
				print "\n";
				if($self->output_file)
				{
					print OUTFILE "Sample\t$stats_header";
					print OUTFILE "\t$alignsum_header" if($self->sample_alignsum_files);
					print OUTFILE "\n";
				}

				$header_printed = 1;
			}

			my ($cov_20x) = split(/\t/, $sample_stats);

			print "$sample_name\t$sample_stats";
			print "\t" . $alignment_summaries{$sample_name} if($alignment_summaries{$sample_name});
			print "\tFLAG" if($cov_20x < 70);
			print "\n";
			if($self->output_file)
			{
				print OUTFILE "$sample_name\t$sample_stats";
				print OUTFILE "\t" . $alignment_summaries{$sample_name} if($alignment_summaries{$sample_name});
				print OUTFILE "\n" 				
			}

		}


	}

	close($input);

	close(OUTFILE) if($self->output_file);
	
	print $stats{'num_samples'} . " samples in file\n";
	print $stats{'with_stats'} . " have RefCov reports\n";
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




#############################################################
# ParseFile - takes input file and parses it
#
#############################################################

sub parse_stats_file
{
	(my $stats_file) = @_;
	
	## Parse the file ##

	my $input = new FileHandle ($stats_file);
	my $lineCounter = 0;
	
	my %coverage_by_depth = ();
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		if($lineCounter > 1)
		{
			my @lineContents = split(/\t/, $line);
			my $num_targets = $lineContents[0];
			my $depth = $lineContents[1];
			my $pct_targets_cov_80 = $lineContents[13];

			$coverage_by_depth{$depth} = $pct_targets_cov_80;
		}
		
	}
	
	close($input);

	my $coverage_header = my $coverage_result = "";

	foreach my $depth (sort descending keys %coverage_by_depth)
	{
		$coverage_header .= "\t" if($coverage_header);
		$coverage_header .= $depth . "x";
		
		$coverage_result .= "\t" if($coverage_result);
		$coverage_result .= $coverage_by_depth{$depth};
	}
	
	sub descending
	{
		$b <=> $a;
	}

	my $coverage_report = $coverage_header . "\n" . $coverage_result;

	return($coverage_report);
}



#############################################################
# ParseFile - takes input file and parses it
#
#############################################################

sub load_alignment_summaries
{
	my ($FileName) = shift(@_);

	my %sample_results = ();

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	my $header_printed = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		(my $sample_name, my $alignsum_file) = split(/\t/, $line);

		if(-e $alignsum_file)
		{
			(my $alignment_summary) = parse_alignsum_file($alignsum_file);
			$sample_results{$sample_name} = $alignment_summary;
		}


	}

	close($input);
	
	return(%sample_results);
}



#############################################################
# ParseFile - takes input file and parses it
#
#############################################################

sub parse_alignsum_file
{
	(my $stats_file) = @_;
	
	## Parse the file ##

	my $input = new FileHandle ($stats_file);
	my $lineCounter = 0;

	my @headers = ();
	
	my %summary = ();
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		my @lineContents = split(/\t/, $line);
		my $numContents = @lineContents;
		
		if($lineCounter == 1)
		{
			for(my $colCounter = 0; $colCounter < $numContents; $colCounter++)
			{
				$headers[$colCounter] = $lineContents[$colCounter];
			}
		}
		elsif($lineCounter > 1)
		{
			for(my $colCounter = 0; $colCounter < $numContents; $colCounter++)
			{
				my $field_name = $headers[$colCounter];
				my $field_value = $lineContents[$colCounter];
				$summary{$field_name} = $field_value;
			}
		}
		
	}
	
	close($input);

	my $on_target_rate = ($summary{'unique_target_aligned_bp'} + $summary{'duplicate_target_aligned_bp'}) / ($summary{'unique_target_aligned_bp'} + $summary{'duplicate_target_aligned_bp'} + $summary{'unique_off_target_aligned_bp'} + $summary{'duplicate_off_target_aligned_bp'});
	$on_target_rate = sprintf("%.2f", $on_target_rate * 100);

	## Change all keys to gigabases ##
	
	foreach my $key (sort keys %summary)
	{
		my $adjusted_value = $summary{$key} / 1000000000;
		$adjusted_value = sprintf("%.2f", $adjusted_value);
		$summary{$key} = $adjusted_value;
	}

	my $result = join("\t", $summary{'unique_target_aligned_bp'}, $summary{'duplicate_target_aligned_bp'}, $summary{'unique_off_target_aligned_bp'}, $summary{'duplicate_off_target_aligned_bp'}, $on_target_rate);
	
	
	
	
	return($result);
}

1;

