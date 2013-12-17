
package Genome::Model::Tools::Ssaha::AlignToGenome;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# AlignToGenome.pm - 	Align reads to a reference genome using ssaha
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

class Genome::Model::Tools::Ssaha::AlignToGenome {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		query_file	=> { is => 'Text', doc => "454 reads in FASTA format" },
		output_file	=> { is => 'Text', doc => "Output file for Ssaha2 alignments" },
                reference	=> { is => 'Text', doc => "Path to ssaha-indexed reference [Defaults to Hs36]", is_optional => 1 },
                ssaha2_params	=> { is => 'Text', doc => "Options for ssaha2 [-454 -best 1 -udiff 1 -output sam]", is_optional => 1 },
                cdna	=> { is => 'Text', doc => "If set to 1, applies Conrad parameters for cDNA alignment [-454 -seeds 2 -diff -1 -kmer 12 -skip 5]", is_optional => 1 },		
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Align reads to a reference genome using Ssaha"                 
}

sub help_synopsis {
    return <<EOS
This command retrieves the locations of unplaced reads for a given genome model
EXAMPLE:	gmt ssaha --query-file s_1_sequence.fastq --output-file s_1_sequence.Hs36.ssaha
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
	my $query_file = $self->query_file;
	my $output_file = $self->output_file;

	## Set default ssaha2 params ##

	my $ssaha2_params = "-454 -best 1 -udiff 1 -output sam";
	$ssaha2_params = $self->ssaha2_params if($self->ssaha2_params);



	if(!(-e $query_file))
	{
		die "Error: Query file not found!\n";
	}

	## Define Ssaha Reference (default to Hs36)

	my $reference = Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.ssaha2';
#	my $reference = Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.ssaha2';	
	my $reference_fasta = Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa';

        if(defined($self->reference))
	{
#		if(-e $self->reference)
#		{
			$reference = $self->reference;
			$reference_fasta = $self->reference;
#		}
#		else
#		{
#			die "Error: Reference file not found!\n";
#		}
	}

	## Build FASTQ File ##
	my $fq_file = "$query_file.fq";

	if(!(-e $fq_file))
	{
		print "Converting FASTA to FASTQ format...\n";
		system("perl ~dkoboldt/Software/SSAHA/fastq.pl $query_file $query_file");
	
		print "Calculating FASTQ...\n";
		system("/gscuser/dkoboldt/Software/SSAHA/other_codes/x86_64/compFastq $query_file.fastq $query_file.fq");
	}

	print "SSAHA FQ file: $fq_file\n";

	if($self->cdna)
	{
		print "Aligning $fq_file to $reference\n";
		$ssaha2_params = "-454 -seeds 2 -diff -1 -kmer 12 -skip 5 -output sam";
		system("bsub -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT} -R\"select[mem>8000] rusage[mem=8000]\" -M 8000000 -oo $output_file.log ssaha2 $ssaha2_params -outfile $output_file $reference_fasta $fq_file");				
	}
	else
	{
		print "Aligning $fq_file to $reference\n";
		system("bsub -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT} -R\"select[mem>8000] rusage[mem=8000]\" -M 8000000 -oo $output_file.log ssaha2 $ssaha2_params -outfile $output_file $reference_fasta $fq_file");		
	}




	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;

