
package Genome::Model::Tools::Analysis::Xenograft::GetMouseAligned;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# GetMouseAligned - Get reads that aligned to the mouse genome
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	07/14/2009 by D.K.
#	MODIFIED:	07/14/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::Xenograft::GetMouseAligned {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		map_file	=> { is => 'Text', doc => "Map file of mouse alignments", is_optional => 0 },
		output_file	=> { is => 'Text', doc => "Output file for read names", is_optional => 0 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Gets a list of reads that aligned to the mouse genome"                 
}

sub help_synopsis {
    return <<EOS
This command gets a list of reads that aligned to the mouse genome
EXAMPLE:	gmt analysis Xenograft get-mouse-aligned --map-file my.map --output-file my.reads
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
	my $map_file = $self->map_file;
	my $output_file = $self->output_file;

	system("maq mapview $map_file | cut -f 1 >$output_file");
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}





sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}


1;

