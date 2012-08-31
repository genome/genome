
package Genome::Model::Tools::Bowtie::AlignToGenome;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# AlignToGenome.pm - 	Align reads to a reference genome using Bowtie
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
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

## Bowtie Parameters ##
my $batch_size = 1000000;
my $num_cores = 1;
my $lsf_queue = "long";

my $novoalign_params = "-c $num_cores -a -l 36 -t 240 -k";	# -o SAM

my $path_to_novoalign = "/gscuser/dkoboldt/Software/NovoCraft/novocraftV2.05.13/novocraft/novoalign";
my $novoalign_reference = Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.novoindex-k14-s3-v2.05.13';




class Genome::Model::Tools::Bowtie::AlignToGenome {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		query_file	=> { is => 'Text', doc => "Illumina/Solexa reads in FASTQ format" },
		output_file	=> { is => 'Text', doc => "Output file for Bowtie alignments" },
                reference	=> { is => 'Text', doc => "Path to bowtie-indexed reference [Defaults to Hs36]", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Align reads to a reference genome using Bowtie"                 
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

	## Get required parameters ##
	my $query_file = $self->query_file;
	my $output_file = $self->output_file;

	if(!(-e $query_file))
	{
		die "Error: Query file not found!\n";
	}

	## Define Bowtie Reference (default to Hs36)

	my $reference = $novoalign_reference;

        if(defined($self->reference))
	{
		if(-e $self->reference)
		{
			$reference = $self->reference;
		}
		else
		{
			die "Error: Reference file not found!\n";
		}
	}

	system("bsub -q $lsf_queue -R\"select[type==LINUX64 && model != Opteron250 && mem>12000] rusage[mem=12000] span[hosts=1]\" -n $num_cores -M 12000000 \"$path_to_novoalign $novoalign_params -d $reference -f $query_file >$output_file\"");

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;

