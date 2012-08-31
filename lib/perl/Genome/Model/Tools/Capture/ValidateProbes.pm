package Genome::Model::Tools::Capture::ValidateProbes;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ValidateProbes - Validate a Probeset with 454 Data
#					
#	AUTHOR:		Will "KingKong" Schierding (wschierd@genome.wustl.edu)
#
#	CREATED:	1/28/2009 by W.S.
#	MODIFIED:	1/28/2009 by W.S.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Capture::ValidateProbes {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		barcode	=> { is => 'Text', doc => "Barcode (usually found from Wiki)", is_optional => 0 },
		sff_file	=> { is => 'Text', doc => "sff file for the given project", is_optional => 0 },
		output_dir	=> { is => 'Text', doc => "Output Directory for files (example: runmap_out)" , is_optional => 0},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Validate a Probeset with 454 Data"                 
}

sub help_synopsis {
    return <<EOS
Create and build genome models for capture datasets
EXAMPLE: gmt capture validate-probes --barcode [required] --sff-file [required] --output-dir [required]
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
	my $barcode = $self->barcode;
	my $sff_file = $self->sff_file;
	my $output_dir = $self->output_dir;

	my $cmd_line = "capture_file_dumper --barcode $barcode --output-type oligo-fasta --output-file capture-probes.fasta";
	system ($cmd_line);

	my $cmd_line2 = "runMapping -pairt -o $output_dir capture-probes.fasta $sff_file";
	system ($cmd_line2);
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


1;

