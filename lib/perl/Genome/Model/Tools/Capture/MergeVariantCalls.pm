
package Genome::Model::Tools::Capture::MergeVariantCalls;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# MergeVariantCalls - Build Genome Models for Capture Datasets
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	12/09/2009 by D.K.
#	MODIFIED:	12/09/2009 by D.K.
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

class Genome::Model::Tools::Capture::MergeVariantCalls {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		varscan_file	=> { is => 'Text', doc => "File of variants in Varscan format", is_optional => 0, is_input => 1 },
		glf_file	=> { is => 'Text', doc => "File of variants in glfSomatic format", is_optional => 0, is_input => 1 },
		output_file	=> { is => 'Text', doc => "Output file to contain merged results" , is_optional => 0, is_input => 1, is_output => 1},
	],
	
	has_param => [
		lsf_resource => { default_value => 'select[model!=Opteron250 && type==LINUX64 && mem>6000] rusage[mem=6000]'},
       ],	
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Merges Varscan and glfSomatic variant calls"                 
}

sub help_synopsis {
    return <<EOS
Merges Varscan and glfSomatic variant calls
EXAMPLE:	gmt capture merge-variant-calls ...
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

	my $varscan_file = $self->varscan_file;
	my $glf_file = $self->glf_file;
	my $output_file = $self->output_file;
	
	if(!(-e $varscan_file && -e $glf_file))
	{
		die "One or more files didn't exist!\n";
	}

	## Run the merge using Varscan ##
	
	my $cmd = "java -Xms3000m -Xmx3000m -jar /gsc/scripts/lib/java/VarScan/VarScan.jar compare $varscan_file $glf_file merge $output_file.unsorted";
	system($cmd);

	$cmd = "gmt capture sort-by-chr-pos --input-file $output_file.unsorted --output-file $output_file";
	system($cmd);
	
	system("rm -rf $output_file.unsorted") if(-e "$output_file");


	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}






1;

