
package Genome::Model::Tools::BwaSw::AlignToGenome;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# AlignToGenome.pm - 	Align reads to a reference genome using bwa smith waterman alignment
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

class Genome::Model::Tools::BwaSw::AlignToGenome {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		query_file	=> { is => 'Text', doc => "454 reads in FASTA format" },
		output_file	=> { is => 'Text', doc => "Output file for Ssaha2 alignments" },
                reference	=> { is => 'Text', doc => "Path to ssaha-indexed reference [Defaults to Hs36]", is_optional => 1 },
                bwa_params	=> { is => 'Text', doc => "Options for bwa-sw [none]", is_optional => 1 },
                bwa_path	=> { is => 'Text', doc => "Path to BWA version to use", is_optional => 0, default => '/gsc/pkg/bio/bwa/bwa-0.5.9-64/bwa' },
		
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Align reads to a reference genome using bwa-sw"                 
}

sub help_synopsis {
    return <<EOS
This command aligns reads to the reference genome using bwa-sw
EXAMPLE:	gmt bwa-sw --query-file s_1_sequence.fastq --output-file s_1_sequence.Hs36.bwa-sw.sam
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

	my $bwa_params = "";
	$bwa_params = $self->bwa_params if($self->bwa_params);



	if(!(-e $query_file))
	{
		die "Error: Query file not found!\n";
	}

	## Define Ssaha Reference (default to Hs36)

	my $reference = Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa';

        if(defined($self->reference))
	{
#		if(-e $self->reference)
#		{
			$reference = $self->reference;
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

	print "FQ file: $fq_file\n";

	print "Aligning $fq_file to $reference\n";
	system("bsub -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER} -R\"select[mem>8000] rusage[mem=8000]\" -M 8000000 -oo $output_file.log \"" . $self->bwa_path . " bwasw $bwa_params $reference $fq_file > $output_file\"");
#	system("ssaha2 $bwa_params -outfile $output_file -save $reference $fq_file");
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;

