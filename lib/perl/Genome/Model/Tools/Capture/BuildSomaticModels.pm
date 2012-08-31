
package Genome::Model::Tools::Capture::BuildSomaticModels;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# BuildSomaticModels - Compare tumor versus normal models to find somatic events
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

class Genome::Model::Tools::Capture::BuildSomaticModels {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		processing_profile	=> { is => 'Text', doc => "Processing profile to use", is_optional => 0, is_input => 1, default => 'Jan 2012 Default Somatic Variation' },
		previously_discovered_variations_build	=> { is => 'Text', doc => "dbSNP list to use for novel filter [default whitelist-135]", is_optional => 0, is_input => 1, default => 110108854 },	
		annotation_build	=> { is => 'Text', doc => "Annotation build to use [default b37]", is_optional => 0, is_input => 1, default => 106409619 },	
		paired_samples_list	=> { is => 'Text', doc => "Tab-delimited file of tumor sample, normal-model-id, tumor-model-id, [opt normal sample]" , is_optional => 0},
		model_basename	=> { is => 'Text', doc => "String to use for naming models; sample will be appended" , is_optional => 0},
		report_only	=> { is => 'Text', doc => "Flag to skip actual execution" , is_optional => 1},
		use_bsub	=> { is => 'Text', doc => "If set to 1, will submit define command to short queue" , is_optional => 1},
		start_models	=> { is => 'Text', doc => "If set to 1, will start models after finding/creating them" , is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Define and build somatic-variation pipeline models for tumor-normal pairs"                 
}

sub help_synopsis {
    return <<EOS
Define and build somatic-variation pipeline models for tumor-normal pairs
EXAMPLE:	gmt capture build-somatic-models --paired-samples-list Paired-Model-IDs.tsv --model-basename "MySomaticVariation-MyPP"
	Paired Model IDs should be formatted like this (no header, and tumor_sample must be valid subject_name):
	tumor_sample	normal_model_id	tumor_model_id normal_sample	
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
	my $processing_profile = $self->processing_profile;
	my $annotation_build = $self->annotation_build;
	my $dbsnp_build = $self->previously_discovered_variations_build;
	my $sample_list = $self->paired_samples_list;
	my $model_basename = $self->model_basename;

	## Parse the sample file ##

	my $input = new FileHandle ($sample_list);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		(my $tumor_sample_name, my $normal_model_id, my $tumor_model_id, my $normal_sample_name) = split(/\t/, $line);
		$stats{'num_pairs'}++;

		my @temp = split(/\-/, $tumor_sample_name);
		my $patient_id = join("-", $temp[0], $temp[1], $temp[2]);
		
		if($normal_sample_name)
		{
			$normal_sample_name =~ s/\$patient_id\-//;
		}

		my $model_name = $model_basename . "-" . $tumor_sample_name;
		$model_name .= "_" . $normal_sample_name if($normal_sample_name);

		my $model_id = get_model_id($model_name);

		print "$tumor_sample_name\t$model_name\n";

		## Build the somatic model ##
		if(!$model_id)
		{
			my $cmd = "genome model define somatic-variation --processing-profile-name \"$processing_profile\" --previously-discovered-variations-build $dbsnp_build --annotation-build $annotation_build --subject-name \"$tumor_sample_name\" --model-name \"$model_name\" --normal-model $normal_model_id --tumor-model $tumor_model_id";
			if($self->use_bsub)
			{
				system("bsub -q short $cmd") if(!$self->report_only);				
			}
			else
			{
				system("$cmd") if(!$self->report_only);
				$model_id = get_model_id($model_name) if(!$self->report_only);
			}
		}

		if($self->start_models && $model_id)
		{
			my $cmd = "genome model build start $model_id";
			print "RUN: $cmd\n";

			if(!$self->report_only)
			{
				system($cmd);
			}
		}

	}

	close($input);
	
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





1;

