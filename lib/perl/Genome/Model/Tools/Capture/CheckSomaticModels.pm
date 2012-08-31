
package Genome::Model::Tools::Capture::CheckSomaticModels;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# CheckSomaticModels - Compare tumor versus normal models to find somatic events
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	12/09/2009 by D.K.
#	MODIFIED:	08/13/2010 by D.K.
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

class Genome::Model::Tools::Capture::CheckSomaticModels {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		processing_profile	=> { is => 'Text', doc => "Processing profile to use [Somatic-Capture-NoSV-Tier1only-Map40-Score40]", is_optional => 1 },
		data_dir	=> { is => 'Text', doc => "Output directory for comparison files" , is_optional => 0},
		sample_list	=> { is => 'Text', doc => "Text file of sample, normal-model-id, tumor-model-id" , is_optional => 0},
		subject_type	=> { is => 'Text', doc => "Subject type, e.g. sample_name, library_name [library_name]" , is_optional => 1},
		model_basename	=> { is => 'Text', doc => "String to use for naming models; sample will be appended" , is_optional => 0},
		report_only	=> { is => 'Text', doc => "Flag to skip actual execution" , is_optional => 1},
		use_bsub	=> { is => 'Text', doc => "If set to 1, will submit define command to short queue" , is_optional => 1},
		build_mafs	=> { is => 'Text', doc => "If set to 1, will build MAF files for completed models" , is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Check the status of somatic models"                 
}

sub help_synopsis {
    return <<EOS
Check the status of somatic models
EXAMPLE:	gt capture check-somatic-models ...
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
	my $processing_profile = "Somatic-Capture-NoSV-Tier1only-Map40-Score40";
	$processing_profile = $self->processing_profile if($self->processing_profile);

	my $sample_list = $self->sample_list;
	my $subject_type = "library_name";
	$subject_type = $self->subject_type if($self->subject_type);

	my $model_basename = $self->model_basename;

	my $data_dir = "./";
	$data_dir = $self->data_dir if($self->data_dir);

	## Print the header ##
	
	print "NORMAL_SAMPLE_NAME\tTUMOR_SAMPLE_NAME\tMODEL_ID\tBUILD_ID\tSTATUS\tVSCAN\tSNIPER\tMERGED\tFILTER\tNOVEL\tTIER1\tINDELS\tTIER1\n";


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
#		print "Getting model id for $model_name...\n";

		my $model_id = get_model_id($model_name);
#		print "Got model id $model_id\n";
	
		my $model_status = "Unknown";

		## Build the somatic model ##
		if(!$model_id)
		{
			print "$tumor_sample_name\tNo Model Named $model_name\n";
		}
		else
		{
			my $model = Genome::Model->get($model_id);			
			
			my @builds = $model->builds;
			
			foreach my $build (@builds)
			{
				my $build_dir = $build->data_directory;
				my $build_status = $build->status;
				$stats{$build_status}++;

				my %build_results = ();
				## Determine paths to key files ##
							
				$build_results{'varscan_somatic_snv'} = get_file_status("$build_dir/varScan.output.snp.formatted.Somatic.hc");
				$build_results{'sniper_somatic_snv'} = get_file_status("$build_dir/somaticSniper.output.snp.filter.hc.somatic");
				$build_results{'merged_somatic_snv'} = get_file_status("$build_dir/merged.somatic.snp");
				$build_results{'merged_somatic_snv_filt'} = get_file_status("$build_dir/merged.somatic.snp.filter");
				$build_results{'merged_somatic_snv_filt_novel'} = get_file_status("$build_dir/merged.somatic.snp.filter.novel");
				$build_results{'merged_somatic_snv_filt_novel_tier1'} = get_file_status("$build_dir/merged.somatic.snp.filter.novel.tier1");

				print join("\t", $normal_sample_name, $tumor_sample_name, $model_id, $build->id, $build->status) . "\t"; #, $build->data_directory
				print join("\t", $build_results{'varscan_somatic_snv'}, $build_results{'sniper_somatic_snv'}, $build_results{'merged_somatic_snv'}) . "\t";
				print join("\t", $build_results{'merged_somatic_snv_filt'}, $build_results{'merged_somatic_snv_filt_novel'}, $build_results{'merged_somatic_snv_filt_novel_tier1'}) . "\t";
				print "\n";
			}

#			exit(0);
			

		}

	}

	close($input);


	print $stats{'snvs_completed'} . " patients have SNVs completed\n";
	print $stats{'indels_completed'} . " patients have indels completed\n";
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}





#############################################################
# ParseFile - takes input file and parses it
#
#############################################################

sub get_file_status
{
	my $file = shift(@_);
	
	if(-e $file)
	{
		my $status = `cat $file | wc -l`;
		chomp($status);
		return($status);
	}
	else
	{
		return("-");
	}
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


#############################################################
# ParseFile - takes input file and parses it
#
#############################################################

sub get_build_ids
{
	my $dir = shift(@_);
	
	my @build_ids = ();
	my $numBuilds = 0;
	
	my $dir_list = `ls -d $dir/build* 2>/dev/null`;
	chomp($dir_list);
	
	my @dir_lines = split(/\n/, $dir_list);
	foreach my $dir_line (@dir_lines)
	{
		my @lineContents = split(/\//, $dir_line);
		my $numContents = @lineContents;
		
		my $short_dir_name = $lineContents[$numContents - 1];
		
		my $build_id = $short_dir_name;
		$build_id =~ s/build//;
		$build_id =~ s/\@//;
		if($build_id)
		{
			$build_ids[$numBuilds] = $build_id;
			$numBuilds++;
		}
	}
	
	return(@build_ids);
}


1;

