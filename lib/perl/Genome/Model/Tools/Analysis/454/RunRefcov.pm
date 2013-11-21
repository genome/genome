
package Genome::Model::Tools::Analysis::454::RunRefcov;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# SamToBam - Align reads with SSAHA2 or other aligner
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	02/25/2009 by D.K.
#	MODIFIED:	02/25/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

my $ref_seq = Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa';
my $ref_index = Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa.fai';

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::454::RunRefcov {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		samples_file	=> { is => 'Text', doc => "Tab-delimited file of sample and SFF file(s)" },
		output_dir	=> { is => 'Text', doc => "Output directory" },
		aligner		=> { is => 'Text', doc => "Aligner that was used." },
		bed_file	=> { is => 'Text', doc => "BED file of RefCov targets." },
		min_depth	=> { is => 'Text', doc => "Min depth filter if other than 1,5,10,15,20", is_optional => 1, default => "1,5,10,15,20" },
		wingspan	=> { is => 'Text', doc => "Wingspan for reporting coverage", is_optional => 1, default => "0" },
		reference		=> { is => 'Text', doc => "Reference sequence [default=Hs36 ssaha2]", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Runs RefCov on target regions using alignment BAM files"                 
}

sub help_synopsis {
    return <<EOS
This command converts SAM files to BAM files
EXAMPLE:	gmt analysis 454 run-refcov --samples-file data/samples.tsv --output-dir data --aligner ssaha2 --bed-file myTargets.bed
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
	my $samples_file = $self->samples_file;
	my $output_dir = $self->output_dir;
	my $bed_file = $self->bed_file;
	my $aligner = $self->aligner;
	my $min_depth = $self->min_depth;
	my $wingspan = $self->wingspan;


	if($self->reference)
	{
		$ref_seq = $self->reference;
		$ref_index = $ref_seq . ".fai";

		die "Reference files not found!\n$ref_seq\n$ref_index\n" if(!(-e $ref_seq && -e $ref_index));
	}

	if(!(-e $samples_file))
	{
		die "Error: Samples file $samples_file not found!\n";
	}

	if(!(-e $bed_file))
	{
		die "Error: BED file $bed_file not found!\n";
	}

	## Parse out bed file shortname ##
	
	my @tempArray = split(/\//, $bed_file);
	my $numContents = @tempArray;
	my $bed_file_short_name = $tempArray[$numContents - 1];


	my $input = new FileHandle ($samples_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		my @lineContents = split(/\t/, $line);			
		my $sample_name = $lineContents[0];
		my $sff_files = $lineContents[1];				
		
		## Create sample output directories ##
		my $sample_output_dir = $output_dir . "/" . $sample_name;
		my $aligner_output_dir = $sample_output_dir . "/" . $aligner . "_out";
		my $aligner_scripts_dir = $sample_output_dir . "/" . $aligner . "_out/scripts";		
		my $refcov_output_dir = $sample_output_dir . "/reference_coverage";
		my $bam_file = "$aligner_output_dir/$sample_name.$aligner.bam";
		
		if(-e $bam_file)
		{
			mkdir($refcov_output_dir) if(!(-d $refcov_output_dir));
			mkdir($aligner_scripts_dir) if(!(-d $aligner_scripts_dir));
			
			## open outfile ##
			my $script_filename = "$aligner_scripts_dir/run_refcov.sh";
			open(SCRIPT, ">$script_filename") or die "Can't open script file: $!\n";
			
			## Step 1: Convert SAM to BAM ##
			print SCRIPT "Running RefCov...\n";
			my $cmd = "gmt bio-samtools ref-cov ";
			$cmd .= "--bed-file $bed_file ";
			$cmd .= "--bam-file $bam_file ";
			$cmd .= "--min-depth-filter 1,5,10,15,20 ";
			$cmd .= "--wingspan $wingspan ";
			$cmd .= "--output-directory $refcov_output_dir ";
			
			print SCRIPT "$cmd\n";


			## Step 2: Generate stats summary ##

			my $stats_file = $refcov_output_dir . "/wingspan_0/$sample_name.$aligner" . "_" . $bed_file_short_name . "_STATS.tsv";
			my $summary_file = $refcov_output_dir . "/wingspan_0/$sample_name.$aligner" . "_" . $bed_file_short_name . "_STATS.txt";
			$cmd = "gmt bio-samtools stats-summary --stats-file $stats_file --output-file $summary_file";
			print SCRIPT "$cmd\n";			

			close(SCRIPT);
			
			## Change permissions ##
			system("chmod 755 $script_filename");
			
			## Run bsub ##
#			system("$script_filename");			
			system("bsub -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER} -R\"select[type==LINUX64 && model != Opteron250 && mem>4000 && tmp>20000] rusage[mem=4000]\" -oo $script_filename.out $script_filename");

		}
		else
		{
			die "Alignment BAM file $bam_file not found!\n";
		}
		
	}
	
	close($input);

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;

