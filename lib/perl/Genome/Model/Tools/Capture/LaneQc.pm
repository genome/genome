
package Genome::Model::Tools::Capture::LaneQc;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# LaneQc - Run per-lane QC for samples in a model group
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	12/23/2010 by D.K.
#	MODIFIED:	12/23/2010 by D.K.
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

class Genome::Model::Tools::Capture::LaneQc {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		group_id		=> { is => 'Text', doc => "ID of model group" , is_optional => 0},
		genotype_files	=> { is => 'Text', doc => "Tab-delimited list of samples and paths to array genotype data", is_optional => 1, is_input => 1 },	
		output_file     => { is => 'Text', doc => "Output file to receive QC results", is_optional => 1, is_input => 1, is_output => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Run per-lane sample QC on succeeded builds in a model group"                 
}

sub help_synopsis {
    return <<EOS
This command runs per-lane sample QC on succeeded builds in a model group
EXAMPLE:	gmt capture lane-qc --group-id 661 --genotype-files Sample-Genotype-Files.tsv
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

	my $group_id = $self->group_id;

	## Open Optional Output Files ##

	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";
	}

	## Keep stats in a single hash ##
	
	my %stats = ();
	
	## Save model ids by subject name ##
	
	my %succeeded_models_by_sample = ();

	## Get the models in each model group ##

	my $model_group = Genome::ModelGroup->get($group_id);
	my @models = $model_group->models; 

	foreach my $model (@models)
	{
		$stats{'models_in_group'}++;
		
		my $model_id = $model->genome_model_id;
		my $subject_name = $model->subject_name;
		
		my $build_dir = my $bam_file = my $snp_file = "";
		my $model_status = "New";


		my $num_builds = 0;		

		my $build_ids = my $build_statuses = "";
		my @builds = $model->builds;

		if($model->last_succeeded_build_directory)
		{
			foreach my $build (@builds)
			{
				my $build_id = $build->id;
				
				print join("\t", $subject_name, $model_id, $build_id, $build->status) . "\n";	
				
				if($build->status eq "Succeeded")
				{
					## Get the BAM File(s)
					
					my @bam_list = $build->get_alignment_bams;
					
					foreach my $bam_file (@bam_list)
					{
						my $instrument_data_id = 0;
						## Parse out instrument data ID ##
						my @pathContents = split(/\//, $bam_file);
						my $numContents = @pathContents;
						
						for(my $colCounter = 0; $colCounter < $numContents; $colCounter++)
						{
							$instrument_data_id = $pathContents[$colCounter + 1] if($pathContents[$colCounter] && $pathContents[$colCounter] eq "alignment_data");
						}
						print "$instrument_data_id\t$bam_file\n";
					}

					
					## Get the SNP file ##
					my $snp_file;

					my $search_string = "ls " . $model->last_succeeded_build_directory . "/snp*/snps_all_sequences.filtered 2>/dev/null | tail -1";
					my $snp_list_result = `$search_string`;
					chomp($snp_list_result) if($snp_list_result);
					
					## Backward compatibility for old names of SNP files ##
					if(!$snp_list_result)
					{
						$search_string = "ls " . $model->last_succeeded_build_directory . "/sam*/filtered.indelpe.snps 2>/dev/null | tail -1";
						$snp_list_result = `$search_string`;
						chomp($snp_list_result) if($snp_list_result);						
					}

					if($snp_list_result && -e $snp_list_result)
					{
						$snp_file = $snp_list_result;
						print "\t$snp_file\n";
						exit(0);
					}

				}
			}
		}

	}	
	
	close(BAMLIST) if($self->output_bam_files);
	close(SNPLIST) if($self->output_snp_files);	

	print $stats{'models_in_group'} . " models in group\n" if($stats{'models_in_group'});
	print $stats{'models_running'} . " models running\n" if($stats{'models_running'});
	print $stats{'models_finished'} . " models finished\n" if($stats{'models_finished'});



	## Determine normal-tumor pairing and completed models ##
	if($self->output_model_pairs)
	{
		my %tumor_sample_names = my %tumor_model_ids = my %normal_model_ids = ();
	
		foreach my $sample_name (keys %succeeded_models_by_sample)
		{
			my $model_id = $succeeded_models_by_sample{$sample_name};
			## Determine patient ID ##
			
			my @tempArray = split(/\-/, $sample_name);
			my $patient_id = join("-", $tempArray[0], $tempArray[1], $tempArray[2]);
			
			## Determine if this sample is normal or tumor ##
			
			my $sample_type = "tumor";
			$sample_type = "normal" if(substr($tempArray[3], 0, 1) eq "1");
	
			if($sample_type eq "tumor")
			{
				$tumor_model_ids{$patient_id} = $model_id;
				$tumor_sample_names{$patient_id} = $sample_name;
			}
			elsif($sample_type eq "normal")
			{
				$normal_model_ids{$patient_id} = $model_id;
			}
			
	#		print "$sample_name\t$patient_id\t$sample_type\n";
		}
		
		foreach my $patient_id (sort keys %tumor_model_ids)
		{
			$stats{'num_patients'}++;
			
			if($normal_model_ids{$patient_id})
			{
				$stats{'num_completed_patients'}++;
	
				my $tumor_sample_name = $tumor_sample_names{$patient_id};
				my $tumor_model_id = $tumor_model_ids{$patient_id};
				my $normal_model_id = $normal_model_ids{$patient_id};
	
				if($self->output_model_pairs)
				{
					print MODELPAIRS "$tumor_sample_name\t$normal_model_id\t$tumor_model_id\n";
				}
	
	#			print "$patient_id\t$tumor_sample_name\t$normal_model_id\t$tumor_model_id\n";
			}
	
		}
		
		close(MODELPAIRS) if($self->output_model_pairs);

		print $stats{'num_patients'} . " patients with models in group\n" if($stats{'num_patients'});
		print $stats{'num_completed_patients'} . " patients with completed tumor+normal builds\n" if($stats{'num_completed_patients'});
		
	}


}




1;

