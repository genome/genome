
package Genome::Model::Tools::Bowtie::BatchCompile;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# BatchCompile.pm - 	Align reads to a reference genome using Bowtie
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	04/22/2009 by D.K.
#	MODIFIED:	04/22/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

## Bowtie Parameters ##
my $batch_size = 1000000;
my $num_cores = 1;
my $lsf_queue = $ENV{GENOME_LSF_QUEUE_BUILD_WORKER};

my %stats = ();

## Bowtie Parameters ##
my $novoalign_params = "-c $num_cores -a -l 36 -t 240 -k";	# -o SAM

my $path_to_novoalign = "/gscuser/dkoboldt/Software/NovoCraft/novocraftV2.05.13/novocraft/novoalign";
my $novoalign_reference = Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.novoindex-k14-s3-v2.05.13';
#my $novoalign_reference = Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.novoindex-k15-s2-v2.05.13';

#my $batch_dir = "/tmp/novoalign";

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Bowtie::BatchCompile {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		batch_id	=> { is => 'Text', doc => "Temporary ID used for naming files" },
		batch_dir	=> { is => 'Text', doc => "Batch directory containing intermediate files"},
		output_file	=> { is => 'Text', doc => "Desired output filename"},
		report_only	=> { is => 'Text', doc => "Report but do not compile results", is_optional => 1},
		resubmit_failures	=> { is => 'Text', doc => "Resubmit failed novoaligns", is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Batch-align reads to a reference genome using Bowtie"                 
}

sub help_synopsis {
    return <<EOS
This command retrieves the locations of unplaced reads for a given genome model
EXAMPLE:	gmt bowtie --query-file s_1_sequence.fastq --output-file s_1_sequence.Hs36.bowtie
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

	
	$stats{'pairs'} = $stats{'pairs_aligned'} = $stats{'reads'} = $stats{'reads_aligned'} = $stats{'unique'} = $stats{'gapped'} = $stats{'qualfilt'} = $stats{'homofilt'} = $stats{'elapsed'} = 0;

	## Get required parameters ##
	my $output_file = $self->output_file;
	my $batch_id = $self->batch_id;
	my $batch_dir = $self->batch_dir;

	if(-e $output_file && !$self->report_only)
	{
		warn "Output file $output_file exists; deleting in 5 seconds\n";
		sleep(5);
		system("rm -rf $output_file");
	}

	## Open the log file ##
	
	open(LOG, ">$output_file.log") or die "Can't open log file: $!\n";
	print LOG "batch\tstatus\tpairs\taligned\treads\taligned\tunique\tgapped\tqualfilt\thomofilt\ttime_elapsed(s)\n";
	print "batch\tstatus\tpairs\taligned\treads\taligned\tunique\tgapped\tqualfilt\thomofilt\ttime_elapsed(s)\n";

	## Find all batch files in the output directory ##
	
	my $batch_files = `ls $batch_dir/$batch_id.*novoalign`;
	chomp($batch_files);
	my @batch_files = split(/\n/, $batch_files);
	foreach my $batch_file (@batch_files)
	{
		my @fileContents = split(/\./, $batch_file);
		my $numContents = @fileContents;
		my $batch_name = $fileContents[$numContents - 2];
		
		$stats{'num_in_batch'}++;
		
		if(has_completed($batch_file))
		{
			my $alignment_stats = alignment_stats($batch_file);
			print LOG $batch_name . "\t" . "Done\t" . $alignment_stats . "\n";
			print $batch_name . "\t" . "Done\t" . $alignment_stats . "\n";
			system("cat $batch_file >>$output_file") if(!$self->report_only);
			$stats{'num_completed'}++;
		}
		elsif(has_failed($batch_file))
		{
			print LOG $batch_name . "\t" . "Failed\n";
			print $batch_file . "\t" . "Failed\n";
			
			$stats{'num_failed'}++;
			
			if($self->resubmit_failures)
			{
				## Infer the FASTQ file ##
				my $fastq_file = $batch_file;
				$fastq_file =~ s/\.novoalign//;
				## Resubmit the failure ##
				system("bsub -q $lsf_queue -R\"select[type==LINUX64 && model != Opteron250 && mem>10000] rusage[mem=10000] span[hosts=1]\" -n $num_cores -M 12000000 -oo $batch_file.log \"$path_to_novoalign $novoalign_params -d $novoalign_reference -f $fastq_file >$batch_file 2>$batch_file.err\"");
			}
		}
		else
		{
			print LOG $batch_name . "\t" . "Running\n";
			print $batch_file . "\t" . "Running\n";
			$stats{'num_running'}++;
		}


	}

	## Print summary stats ##
	
	$stats{'num_in_batch'} = 0 if(!$stats{'num_in_batch'});
	$stats{'num_completed'} = 0 if(!$stats{'num_completed'});
	$stats{'num_failed'} = 0 if(!$stats{'num_failed'});
	$stats{'num_running'} = 0 if(!$stats{'num_running'});

	print "ALL\t";
	print $stats{'num_completed'} . "/" . $stats{'num_in_batch'} . "\t";
	print commify($stats{'pairs'}) . "\t" . commify($stats{'pairs_aligned'}) . "\t" . commify($stats{'reads'}) . "\t" . commify($stats{'reads_aligned'}) . "\t";
	print commify($stats{'unique'}) . "\t" . commify($stats{'gapped'}) . "\t" . commify($stats{'qualfilt'}) . "\t" . commify($stats{'homofilt'}) . "\t" . commify($stats{'elapsed'}) . "\n";

	print $stats{'num_in_batch'} . " alignments in batch\n";
	print $stats{'num_running'} . " are running\n";
	print $stats{'num_failed'} . " have failed\n";
	print $stats{'num_completed'} . " are completed\n";

	

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub has_completed
{
	my $filename = shift(@_);
	my $tail_done = `tail $filename | grep -c Done`;
	chomp($tail_done);
	
	return($tail_done);
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub has_failed
{
	my $filename = shift(@_);
	my $tail_exit = `cat $filename.log | grep -c Exited`;
	chomp($tail_exit);
	
	return($tail_exit);
}



################################################################################################
# Execute - the main program logic
#
################################################################################################

sub alignment_stats
{
	my $filename = shift(@_);

	my $pairs = my $pairs_aligned = my $reads = my $aligned = my $unique = my $gapped = my $qualfilt = my $homofilt =  my $elapsed = "";	

	my $tail = `tail -100 $filename`;
	chomp($tail);

	my @lines = split(/\n/, $tail);
	
	foreach my $line (@lines)
	{
		my @lineContents = split(/\s+/, $line);
		
		if($line =~ "Paired Reads")
		{
			$pairs = $lineContents[3];
		}
		elsif($line =~ "Pairs Aligned")
		{
			$pairs_aligned = $lineContents[3];
		}		
		elsif($line =~ "Read Sequences")
		{
			$reads = $lineContents[3];
		}		
		elsif($line =~ "Aligned")
		{
			$aligned = $lineContents[2];
		}
		elsif($line =~ "Unique Alignment")
		{
			$unique = $lineContents[3];
		}		
		elsif($line =~ "Gapped Alignment")
		{
			$gapped = $lineContents[3];
		}		
		elsif($line =~ "Quality Filter")
		{
			$qualfilt = $lineContents[3];
		}
		elsif($line =~ "Homopolymer Filter")
		{
			$homofilt = $lineContents[3];
		}				
		elsif($line =~ "Elapsed Time")
		{
			($elapsed) = split(/\,/, $lineContents[3]);
		}		
	}

	my $stats_line = "";
	if($pairs)
	{
		$stats_line = commify($pairs) . "\t" . commify($pairs_aligned) . "\t" . commify($reads) . "\t" . commify($aligned) . "\t" . commify($unique) . "\t" . commify($gapped) . "\t" . commify($qualfilt) . "\t" . commify($homofilt) . "\t" . commify($elapsed);		
	}
	else
	{
		$stats_line = "-\t-\t" . commify($reads) . "\t" . commify($aligned) . "\t" . commify($unique) . "\t" . commify($gapped) . "\t" . commify($qualfilt) . "\t" . commify($homofilt) . "\t" . commify($elapsed);
	}

	## Update overall stats ##
	
	#$stats{'pairs'} = $stats{'pairs_aligned'} = $stats{'reads'} = $stats{'reads_aligned'} = $stats{'unique'} = $stats{'gapped'} = $stats{'qualfilt'} = $stats{'elapsed'} = 0;
	$stats{'pairs'} += $pairs;
	$stats{'pairs_aligned'} += $pairs_aligned;
	$stats{'reads'} += $reads;
	$stats{'reads_aligned'} += $aligned;
	$stats{'unique'} += $unique;
	$stats{'gapped'} += $gapped;
	$stats{'qualfilt'} += $qualfilt;
	$stats{'homofilt'} += $homofilt;
	$stats{'elapsed'} += $elapsed;
	
	return($stats_line);
}


sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}

1;

