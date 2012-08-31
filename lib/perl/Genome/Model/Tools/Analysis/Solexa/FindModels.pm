
package Genome::Model::Tools::Analysis::Solexa::FindModels;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# FindModels - Search for specific genome models and check on their status
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	07/24/2009 by D.K.
#	MODIFIED:	07/24/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::Solexa::FindModels {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		sample_list	=> { is => 'Text', doc => "Accept a list of sample names as input", is_optional => 1 },
		sample_name	=> { is => 'Text', doc => "Search by sample name", is_optional => 1 },
		processing_profile_id	=> { is => 'Text', doc => "Search by processing-profile-id", is_optional => 1 },
		print_location	=> { is => 'Text', doc => "If set to 1, prints data location" , is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Searches for existing genome models and checks their status"                 
}

sub help_synopsis {
    return <<EOS
This command searches for Illumina/Solexa data using the database
EXAMPLE:	gmt analysis solexa search-runs --flowcell_id 302RT
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
	my $sample_name = $self->sample_name;
	my $sample_list = $self->sample_list;
	my $processing_profile_id = $self->processing_profile_id;

	my $print_location;
	$print_location = $self->print_location if($self->print_location);


	my $model_list = "";
	my $filter = "";

	## Build the filter based on user input ##

	if($sample_name)
	{
		$filter .= "," if($filter);
		$filter .= "subject_name=$sample_name";
	}

	if($processing_profile_id)
	{
		$filter .= "," if($filter);
		$filter .= "processing_profile_id=$processing_profile_id";		
	}

	$model_list = `genome model list --filter=$filter`;
	chomp($model_list);

	my @lines = split(/\n/, $model_list);
	
	foreach my $line (@lines)
	{
		print "$line\n";
	}

	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}


1;

