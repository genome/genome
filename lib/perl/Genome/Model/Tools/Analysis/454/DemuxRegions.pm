
package Genome::Model::Tools::Analysis::454::DemuxRegions;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# DemuxRegions - Load 454 reads from a sample-SFF tab-delimited file
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

class Genome::Model::Tools::Analysis::454::DemuxRegions {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		sff_list	=> { is => 'Text', doc => "A file listing full paths to the main (non-demuxed) 454 SFF data files" },
		output_file	=> { is => 'Text', doc => "Output file of sample-library-primer-sff" },
		skip_if_present => { is => 'Text', doc => "Skip if SFF/Fasta/Qual files are present", is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "De-multiplexes regions by querying the database"                 
}

sub help_synopsis {
    return <<EOS
This command links de-multiplexed regions to sample names by querying the database.
The file provided to --sff-list should contain paths to one SFF data file per run. This is
usually provided in the e-mail from production and should NOT point to a file in the "demux"
subdirectory. For example:
/gscmnt/sata890/production/105953462/R_2010_12_08_14_36_53_FLX10060105_adminrig_105953462/D_2010_12_09_04_22_29_blade8-4-5_fullProcessing/sff/GSG1JRV01.sff
It is only used to get the path to the SFF data directory, which is then searched for de-multiplexed
files under the "demux" subdirectory. 
EXAMPLE:	gmt analysis 454 demux-regions --sff-list data/run-regions.fof
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
	my $sff_list = $self->sff_list;
	my $output_file = $self->output_file;

	if(!(-e $sff_list))
	{
		die "Error: Samples file not found!\n";
	}

	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	print OUTFILE "run_name\tregion_id\tprimer_seq\tlibrary_name\tsample_name\tsff_file\n";

	## Parse the SFF file list ##

	my $input = new FileHandle ($sff_list);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		if($line && substr($line, length($line) - 4, 4) eq ".sff")
		{
			## Parse out region number ##
			my $region = substr($line, length($line) - 5, 1);

			## Parse out run name ##
			my @lineContents = split(/\//, $line);
			my $run_name = "";
			my $run_directory = "";
			my $region_name = "";

			foreach my $string (@lineContents)
			{
				if($string && substr($string, 0, 2) eq "R_")
				{
					$run_name = $string;
				}
				
				if($string && !($string eq "sff" || $string =~'\.sff'))
				{
					$run_directory .= "/" . $string;
				}
				
				$region_name = $string; ## Last string should be region name ##
			}

			$run_directory =~ s/\s+//g;
			$region_name =~ s/\.sff//;

			## Die if no run name or region parsed is unexpected ##

			if(!$run_name)
			{
				die "Unable to parse a run name from $line\n";
			}
			elsif(!$region_name)
			{
				die "No region name parsed from $line\n";
			}
			elsif($region ne "1" && $region ne "2" && $region ne "3" && $region ne "4")
			{
				die "Unrecognized region: $region from $line\n" ;				
			}
			
			
			## Verify presence of demux_dir ##
			
			my $demux_dir = $run_directory . "/sff/demux";
			
			if(!(-d $demux_dir))
			{
				warn "WARNING: No demux dir found in $run_directory/sff, skipping...\n";
			}
			elsif($run_name && $region)
			{
				warn "Processing run $run_name region $region...\n";
				my $demux_result = process_region($run_name, $region);
				
				my @results = split(/\n/, $demux_result);
				foreach my $result_line (@results)
				{
					my ($primer, my $library_name) = split(/\t/, $result_line);

					## Get the sample name ##
					my $sample_name = get_sample_name($library_name);
					
					## Get the SFF file ##
					
					my $sff_file_path = $demux_dir . "/" . $region_name . ".demux." . $primer . ".sff";

					if(-e $sff_file_path)
					{
						print "$sample_name\t$library_name\t$primer\t$sff_file_path\n";
						print OUTFILE "$run_name\t$region\t$primer\t$library_name\t$sample_name\t$sff_file_path\n";						
					}
					else
					{
						warn "WARNING: $sff_file_path not found!\n";
					}

				}

			}
			else
			{
				warn "WARNING: Unable to parse run name ($run_name) or region ($region) from $line\n";
			}
		}
	}
	
	close($input);

	close(OUTFILE);

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



################################################################################################
# Execute - the main program logic
#
################################################################################################

sub process_region
{
	my ($run_name, $region) = @_;

	my $demux_result = "";
	
	print "$run_name\t$region\n";
	
	my $cmd = "454info --run-name \"$run_name\" --report \"region_index\"";
	my $index_info = `$cmd`;

	my @index_lines = split(/\n/, $index_info);
	foreach my $line (@index_lines)
	{
		my @lineContents = split(/\s+/, $line);
		my $numContents = @lineContents;
		
		if($numContents >= 5)
		{
			## Make sure the region number matches ##
			if($lineContents[1] && $lineContents[1] eq $region)
			{
				## Make sure it's a primer sequence ##
				if($lineContents[2] && length($lineContents[2]) >= 10 && length($lineContents[2]) <= 13)
				{
					my $primer_sequence = $lineContents[2];
					my $library_name = $lineContents[3];
					
					if($library_name)
					{
						$demux_result .= "\n" if($demux_result);
						$demux_result .= "$primer_sequence\t$library_name";
					}
				}
			}
		}
	}

	print "$demux_result\n";
	
	if(!$demux_result)
	{
		warn "WARNING: No demux possible for $run_name\t$region\n";
		warn "454info command resulted in: $index_info\n";
	}
	
	return($demux_result);
}



################################################################################################
# Execute - the main program logic
#
################################################################################################

sub get_sample_name
{
	my $library_name = shift(@_);
	my $cmd = "genome library list --filter=name='$library_name' --noheaders --show=sample_name | head -1";
	my $sample_name = `$cmd`; 
	chomp($sample_name);
	return($sample_name);
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub parse_file
{
	my $FileName = shift(@_);
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	


	}
	
	close($input);
	
}


1;

