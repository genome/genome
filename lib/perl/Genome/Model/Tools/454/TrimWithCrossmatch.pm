
package Genome::Model::Tools::454::TrimWithCrossmatch;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# TrimWithCrossmatch.pm -	Trim reads for M13/MID primer sequence using cross_match alignments, which are more sensitive
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	10/28/2008 by D.K.
#	MODIFIED:	10/28/2008 by D.K.
#
#	NOTES:		
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::454::TrimWithCrossmatch {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		in_sff_file	=> { is => 'Text', doc => "Input SFF file to be trimmed" },
		out_sff_file	=> { is => 'Text', doc => "Output trimmed SFF file to be created" },		
		primer_fasta	=> { is => 'Text', doc => "Fasta file of primer/adaptor sequence to trim", is_optional => 1 },	
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Trim SFF files for vector/primer sequence"                 
}

sub help_synopsis {
    return <<EOS
This command uses cross_match to identify vector/primer sequence at the *ends* (first/last 20 bp) of 454 reads.
The output SFF file contains updated trimming information.
EXAMPLE 1:	

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
	my $in_sff_file = $self->in_sff_file;
	my $out_sff_file = $self->out_sff_file;	
	
	## Set defaults for and then get optional parameters ##
	
	my $primer_fasta = "454-M13-MID.fasta";
	
	$primer_fasta = $self->primer_fasta if($self->primer_fasta);	


	## Verify that input file exists ##
	
	if(!(-e $in_sff_file))
	{
		print "Input SFF file does not exist. Exiting...\n";
		return(0);
	}	

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}





1;

