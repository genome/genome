
package Genome::Model::Tools::Capture::ReportCoverage;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ReportCoverage - Compare tumor versus normal models to find somatic events
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

class Genome::Model::Tools::Capture::ReportCoverage {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		sample_list	=> { is => 'Text', doc => "Text file normal-tumor sample pairs to include, one pair per line" , is_optional => 0},
		coverage_dir	=> { is => 'Text', doc => "Output directory for the ref-cov reports" , is_optional => 0},
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
	my $sample_list = $self->sample_list;
	my $coverage_dir = $self->coverage_dir;

	## Reset statistics ##
	$stats{'num_samples'} = 0;

	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open output file: $!\n";
		print OUTFILE "sample_name\ttotal_regions\tcov_20x\tcov_10x\tcov_8x\tcov_1x\tcov_0x\n";
	}

	print "sample_name\ttotal_regions\tcov_20x\tcov_10x\tcov_8x\tcov_1x\tcov_0x\n";

	my $input = new FileHandle ($sample_list);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		(my $sample_name) = split(/\t/, $line);
		$stats{'num_samples'}++;
	
		if(stats_files_exist($coverage_dir, $sample_name))
		{
			$stats{'with_stats'}++;
			
			my $sample_coverage = report_coverage($coverage_dir, $sample_name);
			print "$sample_name\t$sample_coverage\n";
			print OUTFILE "$sample_name\t$sample_coverage\n" if($self->output_file);
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

sub stats_files_exist
{
	(my $output_dir, my $sample_name) = @_;
	
	## Determine STATS file name ##
	my $coverage = 1;
	my $stats_file_root = $output_dir . "/" . $sample_name . ".stats";

	if(-e "$stats_file_root.1x.tsv" && -e "$stats_file_root.8x.tsv" && -e "$stats_file_root.10x.tsv" && -e "$stats_file_root.20x.tsv")
	{
		return 1;
	}
	
	return 0;
}


#############################################################
# ParseFile - takes input file and parses it
#
#############################################################

sub report_coverage
{
	(my $output_dir, my $sample_name) = @_;
	
	## Determine STATS file name ##

	my $stats_file_root = $output_dir . "/" . $sample_name . ".stats";

	my %cov_stats_1x = parse_stats_file("$stats_file_root.1x.tsv");
	my %cov_stats_8x = parse_stats_file("$stats_file_root.8x.tsv");
	my %cov_stats_10x = parse_stats_file("$stats_file_root.10x.tsv");
	my %cov_stats_20x = parse_stats_file("$stats_file_root.20x.tsv");	

	my $covered_0x = my $covered_1x = my $covered_8x = my $covered_10x = my $covered_20x = 0;

	my $num_regions = $cov_stats_1x{'num_regions'};
	$covered_20x = $cov_stats_20x{'covered_80pct'};
	$covered_10x = $cov_stats_10x{'covered_80pct'} - $covered_20x;
	$covered_8x = $cov_stats_8x{'covered_80pct'} - $covered_20x - $covered_10x;
	$covered_1x = $cov_stats_1x{'covered_80pct'} - $covered_20x - $covered_10x - $covered_8x;
	$covered_0x = $num_regions - $cov_stats_1x{'covered_80pct'};

	my $summary = "$num_regions\t$covered_20x\t$covered_10x\t$covered_8x\t$covered_1x\t$covered_0x";
	return ($summary);
	
	return 0;
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
	
	my %cov_stats = ();
	$cov_stats{'num_regions'} = $cov_stats{'covered_80pct'} = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		(my $region_name, my $pct_covered) = split(/\t/, $_);
		$cov_stats{'num_regions'}++;
		$cov_stats{'covered_80pct'}++ if($pct_covered >= 80);
	}
	
	close($input);

	return(%cov_stats);
}


1;

