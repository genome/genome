
package Genome::Model::Tools::Analysis::Solexa::BuildModels;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# BuildModels - Build models for lanes of Illumina/Solexa data
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	10/16/2009 by D.K.
#	MODIFIED:	10/16/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::Solexa::BuildModels {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		filter	=> { is => 'Text', doc => "Filter for genome instrument-data search", is_optional => 0 },
		pp_name	=> { is => 'Text', doc => "Processing profile name to use for models", is_optional => 0 },
		model_name	=> { is => 'Text', doc => "String to use for naming models [sample will be appended]" , is_optional => 0},
		auto_assign	=> { is => 'Text', doc => "Auto assign instrument data to model", is_optional => 1 },
		auto_build	=> { is => 'Text', doc => "Auto build alignments when data added", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Build models for Illumina lanes"                 
}

sub help_synopsis {
    return <<EOS
This command searches genome instrument-data and builds models for the results
EXAMPLE:	gmt analysis solexa build-models --filter=flow_cell_id=
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
	my $filter = $self->filter;
	my $pp_name = $self->pp_name;
	my $model_name = $self->model_name;


	## Search for instrument data based on filter ##
	
	my $instrument_data = `genome instrument-data list solexa --style=csv --noheaders --filter=$filter --show=id,flow_cell_id,lane,sample_name,library_name`;

	if($instrument_data)
	{
		my @lines = split(/\n/, $instrument_data);
		my $num_lines = @lines;
		
		print "$num_lines Illumina/Solexa lanes returned\n";
		
		foreach my $line (@lines)
		{
			(my $id, my $flow_cell_id, my $lane, my $sample_name, my $library_name) = split(/\,/, $line);
			
			print "ID=$id FC=$flow_cell_id Lane=$lane Sample=$sample_name Lib=$library_name\n";

			## Check to see if model exists ##
			
			my $model_id = get_model_id("$model_name-$sample_name");

			if(!$model_id)
			{
				## Build genome model command ##
				my $cmd = "genome model define reference-alignment ";
				$cmd .= "--processing-profile-name=\"$pp_name\" ";
				$cmd .= "--subject-name=$sample_name ";
				$cmd .= "--model-name=\"$model_name-$sample_name\" ";
				if($self->auto_assign)
				{
					$cmd .= "--auto-assign-inst-data ";
				}
				if($self->auto_build)
				{
					$cmd .= "--auto-build-alignments";
				}
				
				## Execute the define command ##
				
				system($cmd);

				$model_id = get_model_id("$model_name-$sample_name");
			}

			## If model is found, build it ##

			if($model_id)
			{
				print "Model ID:\t$model_id\n";				

				## Assign instrument data to model ##

#				my $cmd = "genome model instrument-data assign expression --model $model_id --instrument-data $id";
				
				## Build the model ##
				
#				$cmd = "genome model build --model-id $model_id";
			}
			


		}
	}
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}





#############################################################
# ParseFile - takes input file and parses it
#
#############################################################

sub get_model_id
{
	my $model_name = shift(@_);
	my $model_id = 0;

	my $model_output = `genome model list --filter=name=\'$model_name\' --show=id 2>/dev/null`;
	chomp($model_output);
	my @output_lines = split(/\n/, $model_output);
	
	foreach my $line (@output_lines)
	{
		$line =~ s/[^0-9]//g;
		if($line)
		{
			$model_id = $line;
		}
	}
	
	return($model_id);
}



sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}


1;

