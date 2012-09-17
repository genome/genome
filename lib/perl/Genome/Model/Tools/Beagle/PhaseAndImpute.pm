package Genome::Model::Tools::Beagle::PhaseAndImpute;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# MutationRate - Calculate the mutation rate (per megabase) given a list of mutations (e.g. tier1 SNVs) and a set of regions (e.g. coding space)
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	04/22/2011 by D.K.
#	MODIFIED:	04/22/2011 by D.K.
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

class Genome::Model::Tools::Beagle::PhaseAndImpute {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		vcf_file	=> { is => 'Text', doc => "Input VCF file (SNVs only)" , is_optional => 0, is_input => '1'},
		output_basename	=> { is => 'Text', doc => "Output basename for BEAGLE files." , is_optional => 1, default => 'beagle', is_input => '1'},
		reference	=> { is => 'Text', doc => "Reference fasta file (defaults to build 37)" , is_optional => 1, default => '/gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa', is_input => '1'},		
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Builds BEAGLE input files and runs phasing/imputation"                 
}

sub help_synopsis {
    return <<EOS
This command builds BEAGLE input files and runs phasing/imputation
EXAMPLE:	gmt beagle phase-and-impute --vcf-file my.vcf --output-file 
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

	my $vcf_file = $self->vcf_file;
	my $reference = $self->reference;
	my $output_basename = $self->output_basename;


	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;
