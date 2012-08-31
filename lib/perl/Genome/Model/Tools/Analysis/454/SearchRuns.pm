
package Genome::Model::Tools::Analysis::454::SearchRuns;     # rename this when you give the module file a different name <--

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

class Genome::Model::Tools::Analysis::454::SearchRuns {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		region_id	=> { is => 'Text', doc => "Search by region_id", is_optional => 1 },
		sample_name	=> { is => 'Text', doc => "Search by sample name", is_optional => 1 },
		run_name	=> { is => 'Text', doc => "Search by run name" , is_optional => 1},
		print_location	=> { is => 'Text', doc => "If set to 1, prints data location" , is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Search for 454 runs using the database"                 
}

sub help_synopsis {
    return <<EOS
This command searches for 454 data using the database
EXAMPLE:	gmt analysis 454 search-runs --region-id 2809629870
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
	my $region_id = $self->region_id;
	my $sample_name = $self->sample_name;
	my $run_name = $self->run_name;
	my $print_location;
	$print_location = $self->print_location if($self->print_location);

	my $sqlrun, my $rows_returned;

	if($region_id)
	{
		$sqlrun = `sqlrun "select region_id, run_name, sample_name, total_key_pass FROM run_region_454 where region_id = $region_id ORDER BY region_id" --instance warehouse --parse`;
	}

	if($sample_name)
	{
#		$sqlrun = `sqlrun "select flow_cell_id, lane, sample_name, library_name, read_length, filt_clusters, seq_id, gerald_directory, median_insert_size, filt_aligned_clusters_pct from 454_lane_summary where sample_name LIKE '\%$sample_name\%' ORDER BY lane" --instance warehouse --parse`;
	}

	if($run_name)
	{
#		$sqlrun = `sqlrun "select flow_cell_id, lane, sample_name, library_name, read_length, filt_clusters, seq_id, gerald_directory, median_insert_size, filt_aligned_clusters_pct from 454_lane_summary where library_name = '$library_name' ORDER BY lane" --instance warehouse --parse`;
	}

	if($sqlrun)
	{
#		print "fcell\tlane\tlibrary_type\tfilt_reads\taln%\tsample_name\tlibrary_name\tstatus\n";
		
		my @lines = split(/\n/, $sqlrun);
		my %lane_pairs = ();
		
		foreach my $line (@lines)
		{
			if($line && (substr($line, 0, 4) eq "FLOW" || substr($line, 0, 1) eq "-"))
			{
				
			}
			elsif($line && $line =~ "Execution")
			{
				($rows_returned) = split(/\s+/, $line);
				print "$rows_returned rows returned\n";
			}
			elsif($line)
			{
				(my $region_id, my $run_name, my $sample_name, my $total_key_pass) = split(/\t/, $line);

				## Get num reads ##
				
#				my $num_reads = commify($filt_clusters);

				
				## Print result ##
				print "$region_id \t$run_name \t$sample_name\t$total_key_pass\n";

			}
		}
	}

	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}


1;

