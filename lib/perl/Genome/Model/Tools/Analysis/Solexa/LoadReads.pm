
package Genome::Model::Tools::Analysis::Solexa::LoadReads;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# LoadReads - Run maq sol2sanger on Illumina/Solexa files in a gerald_directory
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	04/01/2009 by D.K.
#	MODIFIED:	11/03/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;
use Cwd;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::Solexa::LoadReads {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		flowcell_id	=> { is => 'Text', doc => "Search by flowcell_id", is_optional => 1 },
		sample_name	=> { is => 'Text', doc => "Search by sample name", is_optional => 1 },
		library_name	=> { is => 'Text', doc => "Search by library name" , is_optional => 1},
		project_name	=> { is => 'Text', doc => "Search by research project name" , is_optional => 1},
		print_location	=> { is => 'Text', doc => "If set to 1, prints data location" , is_optional => 1},
		include_lanes	=> { is => 'Text', doc => "Specify which lanes of a flowcell to include [e.g. 1,2,3]" , is_optional => 1},
		output_dir	=> { is => 'Text', doc => "Output dir to receive files" , is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Obtains reads from a flowcell_id in FastQ format"                 
}

sub help_synopsis {
    return <<EOS
This command builds FASTQ files for Illumina/Solexa data using the database
EXAMPLE 1:	gmt analysis solexa load-reads --flowcell_id 302RT --include-lanes 1,2,3,4 --output-dir output_dir
EXAMPLE 2:	gmt analysis solexa load-reads --sample-name H_GP-0365n --output-dir H_GP-0365n
EXAMPLE 3:	gmt analysis solexa load-reads --library-name H_GP-0365n-lib2 --output-dir H_GP-0365n
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
	my $flowcell_id = $self->flowcell_id;
	my $sample_name = $self->sample_name;
	my $library_name = $self->library_name;
	my $project_name = $self->project_name;
	my $print_location;
	$print_location = $self->print_location if($self->print_location);
	my $output_dir = "./";
	$output_dir = $self->output_dir if($self->output_dir);

	mkdir($output_dir) if(!(-d $output_dir));

	## Handle include-lanes when specified ##

	my $include_lanes;
	$include_lanes = $self->include_lanes if($self->include_lanes);
	my %lanes_to_include = ();
	
	if($include_lanes)
	{
		my @lanes = split(/\,/, $include_lanes);
		foreach my $desired_lane (@lanes)
		{
			$lanes_to_include{$desired_lane} = 1;
		}
	}
	

	## Get current directory ##
	
	my $cwd = getcwd;

	my $sqlrun, my $rows_returned, my $cmd;

	if($flowcell_id)
	{
		$sqlrun = `sqlrun "select flow_cell_id, lane, sample_name, library_name, read_length, filt_clusters, seq_id, gerald_directory, median_insert_size, filt_aligned_clusters_pct from solexa_lane_summary where flow_cell_id = '$flowcell_id' ORDER BY flow_cell_id, lane" --instance warehouse --parse`;
	}

	if($sample_name)
	{
		$sqlrun = `sqlrun "select flow_cell_id, lane, sample_name, library_name, read_length, filt_clusters, seq_id, gerald_directory, median_insert_size, filt_aligned_clusters_pct from solexa_lane_summary where sample_name LIKE '\%$sample_name\%' ORDER BY lane" --instance warehouse --parse`;
	}

	if($library_name)
	{
		$sqlrun = `sqlrun "select flow_cell_id, lane, sample_name, library_name, read_length, filt_clusters, seq_id, gerald_directory, median_insert_size, filt_aligned_clusters_pct from solexa_lane_summary where library_name = '$library_name' ORDER BY lane" --instance warehouse --parse`;
	}

	if($project_name)
	{
		$sqlrun = `sqlrun "select flow_cell_id, lane, sample_name, library_name, read_length, filt_clusters, seq_id, gerald_directory, median_insert_size, filt_aligned_clusters_pct from solexa_lane_summary where research_project = '$project_name' ORDER BY flow_cell_id, lane" --instance warehouse --parse`;		
	}	

	if($sqlrun)
	{
#		print "$sqlrun\n"; exit(0);
		
		print "fcell\tlane\tlibrary_type\tfilt_reads\taln%\tsample_name\tlibrary_name\tstatus\n";
		
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
				(my $flowcell, my $lane, my $sample, my $library, my $read_length, my $filt_clusters, my $seq_id, my $gerald_dir, my $insert_size, my $align_pct) = split(/\t/, $line);
				
				## Proceed if lane to be included ##
				if(!$include_lanes || $lanes_to_include{$lane})
				{
					## Get num reads ##
					
					my $num_reads = commify($filt_clusters);
					$align_pct = 0 if(!$align_pct);
					$align_pct = sprintf("%.2f", $align_pct) . '%';
					
					## Get SE or PE ##
					
					my $end_type = "SE";
					my $lane_name = $lane;
	
					if($insert_size)
					{
						$end_type = "PE";
						$lane_pairs{"$flowcell.$lane"} = 1 if(!$lane_pairs{"$flowcell.$lane"});
						$lane_name .= "_" . $lane_pairs{"$flowcell.$lane"};
						$lane_pairs{"$flowcell.$lane"}++;
					}
					
					## Create flowcell output dir and fastq output dir if necessary ##
					
					my $flowcell_dir = $output_dir . "/" . $flowcell;
					my $fastq_dir = $output_dir . "/" . $flowcell . "/fastq_dir";
					mkdir($flowcell_dir) if(!(-d $flowcell_dir));
					mkdir($fastq_dir) if(!(-d $fastq_dir));
					my $output_fastq = $fastq_dir . "/" . "s_" . $lane_name . "_sequence.fastq";
					
					## Get gerald dir ##
					
					my $status = "unknown";
					my $location = "";
	
					## FILES STILL ON FILESYSTEM ##
					if(-d $gerald_dir)
					{
						## Try to get the file locations ##
						
						my $find_command = `find $gerald_dir/ -name s_$lane_name\_*sequence.txt`;
						if($find_command)
						{
							chomp($find_command);
							$location = $find_command;
	
							if($location && -e $location)
							{
								my @temp = split(/\//, $location);
								my $sata_drive = $temp[2];
								
								$status = "on_filesystem ($sata_drive)";
								
								## Output the command ##
								
								$cmd = "maq sol2sanger $location $output_fastq";
#								print "$cmd\n";
#								system($cmd);
								system("bsub -q $ENV{GENOME_LSF_QUEUE_SHORT} -R\"select[type==LINUX64 && mem>1000] rusage[mem=1000]\" $cmd");
							}
						}
					}
					
					## EXTRACTING FILES FROM ARCHIVES ##
					if(!($location && -e $location))
					{
						## Try to find the file in the archive ##
						
						my $archive_file = `sqlrun "select path FROM seq_fs_path WHERE seq_id = $seq_id AND data_type = 'illumina fastq tgz'" --instance warehouse | grep sequence_$flowcell`;
						chomp($archive_file) if($archive_file);
						if($archive_file && -e $archive_file)
						{
							$location = $archive_file;
	
							## Move to directory ##
							chdir $flowcell_dir;
							
							## Extract the sequence.txt files ##
							$cmd = "tar xf $archive_file";
							print "$cmd\n";
							my $extracted_files = `$cmd`;
							chomp($extracted_files);
							
							## Extract each file ##
							
							my @files = split(/\n/, $extracted_files);
							foreach my $extracted_file (@files)
							{
								my $this_output_file = $extracted_file;
								$this_output_file =~ s/\.txt/\.fastq/;
								$cmd = "maq sol2sanger $location fastq_dir/$this_output_file";
								print "$cmd\n";

								## Run command, delete .txt if successful ##
								if(!system($cmd))
								{
									system("rm -rf $this_output_file");
								}
							}
							
							## Move back to directory ##
							chdir $cwd;
	
							my @temp = split(/\//, $archive_file);
							my $sata_drive = $temp[2];
							$status = "archived ($sata_drive)";						
						}
					}
					
					## Print result ##
					print "$flowcell \t$lane_name \t$read_length bp $end_type\t$num_reads \t$align_pct \t$sample \t$library \t$status\n";
	
					## Print location ##
					if($print_location && $location)
					{
						print "$location\n";
					}					
				}
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

