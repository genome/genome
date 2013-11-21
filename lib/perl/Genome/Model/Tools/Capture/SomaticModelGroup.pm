
package Genome::Model::Tools::Capture::SomaticModelGroup;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ModelGroup - Build Genome Models for Capture Datasets
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
my %included_variants = ();

my %already_reviewed = ();
my %passed_sites = my %wildtype_sites = my %germline_sites = ();
my %patient_sites_passed_review = ();
my $maf_header = "";
my $maf_header_printed = 0;

class Genome::Model::Tools::Capture::SomaticModelGroup {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		group_id		=> { is => 'Text', doc => "ID of model group" , is_optional => 0},
		output_build_dirs	=> { is => 'Text', doc => "If specified, outputs last succeeded build directory for each sample to this file" , is_optional => 1},
		output_review	=> 	{ is => 'Text', doc => "Specify a directory to output SNV/indel files for manual review" , is_optional => 1},
		output_germline_calls	=> 	{ is => 'Text', doc => "Specify a directory to output germline SNV/indel calls" , is_optional => 1},
		output_loh_calls	=> 	{ is => 'Text', doc => "Specify a directory to output LOH calls" , is_optional => 1},
		process_loh_calls	=> 	{ is => 'Text', doc => "Specify a directory to process LOH calls" , is_optional => 1},
		germline_roi_file	=> 	{ is => 'Text', doc => "A file in BED format to restrict germline calls" , is_optional => 1},
		output_maf_file	=> 	{ is => 'Text', doc => "Output a MAF file for downstream analysis" , is_optional => 1},
		uhc_filter	=> 	{ is => 'Text', doc => "If set to 1, apply ultra-high-conf filter to SNV calls" , is_optional => 1},
		uhc_indels	=> 	{ is => 'Text', doc => "If set to 1, pass hc indels for MAF inclusion" , is_optional => 1},
		reference	=> 	{ is => 'Text', doc => "Reference to use for bam-readcounts-based filters; defaults to build 37" , is_optional => 1, example_values => ['/gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa']},
		reference_transcripts	=> 	{ is => 'Text', doc => "Transcripts to use for variant annotation " , example_values => ['NCBI-human.combined-annotation/58_37c_v2']},
		review_database_snvs	=> 	{ is => 'Text', doc => "If provided, use to exclude already-reviewed sites" , is_optional => 1},
		review_database_indels	=> 	{ is => 'Text', doc => "If provided, use to exclude already-reviewed indels" , is_optional => 1},
		tier_file_location	=> 	{ is => 'Text', doc => "Folder containing tiering BED files (defaults to build 37) " , is_optional => 1, example_values => ['/gscmnt/ams1102/info/model_data/2771411739/build106409619/annotation_data/tiering_bed_files_v3']},
		varscan_copynumber	=> 	{ is => 'Text', doc => "Specify a directory to output VarScan CopyNumber Jobs" , is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Operate on capture somatic model groups"                 
}

sub help_synopsis {
    return <<EOS
Operate on capture somatic model groups
EXAMPLE:	gmt capture somatic-model-group --group-id 3328
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


	## Reset Statistics 
	$stats{'review_snvs_possible'} = $stats{'review_snvs_already'} = $stats{'review_snvs_already_wildtype'} = $stats{'review_snvs_already_germline'} = $stats{'review_snvs_filtered'} = $stats{'review_snvs_included'} = 0;
	$stats{'review_indels_possible'} = $stats{'review_indels_already'} = $stats{'review_indels_already_wildtype'} = $stats{'review_indels_already_germline'} = $stats{'review_indels_filtered'} = $stats{'review_indels_included'} = 0;


	if($self->review_database_snvs)
	{
		print "Loading SNV review...\n";
		load_review_database($self->review_database_snvs);
	}
	
	if($self->review_database_indels)
	{
		print "Loading indel review...\n";
		load_review_database($self->review_database_indels);
	}

	
	## Save model ids by subject name ##
	
	my %succeeded_models_by_sample = ();

	## Output build dirs##
	
	if($self->output_build_dirs)
	{
		open(BUILDDIRS, ">" . $self->output_build_dirs) or die "Can't open outfile: $!\n";
	}

	if($self->output_maf_file)
	{
		open(MAF, ">" . $self->output_maf_file) or die "Can't open MAF file: $!\n";
	}

	## Get the models in each model group ##

	my $model_group = Genome::ModelGroup->get($group_id);
	my @models = $model_group->models;
	my $debug_counter = 0;

	foreach my $model (@models)
	{
		$stats{'models_in_group'}++;
		
		my $model_id = $model->genome_model_id;
		my $model_name = $model->name;
		my $subject_name = $model->subject_name;
		$subject_name = "Model" . $model_id if(!$subject_name);

		## Get normal and tumor model ##
		
		my $normal_model = $model->normal_model;
		my $tumor_model = $model->tumor_model;

		## Get Model IDs ##
		
		my $normal_model_id = $normal_model->id;
		my $tumor_model_id = $tumor_model->id;

		## Get TCGA-Compliant Subject Names ##

		my $normal_subject_name = $normal_model->subject_name;
		my $tumor_subject_name = $tumor_model->subject_name;

		my ($sample_prefix) = split(/\-/, $normal_subject_name);
		$normal_subject_name =~ s/$sample_prefix/TCGA/;
		
		($sample_prefix) = split(/\-/, $tumor_subject_name);
		$tumor_subject_name =~ s/$sample_prefix/TCGA/;

		my @temp = split(/\-/, $tumor_subject_name);
		my $patient_id = join("-", $temp[0], $temp[1], $temp[2]);

		my @normal_builds = $normal_model->builds;
		my @tumor_builds = $tumor_model->builds;

		## Get normal and tumor builds ##
		my $normal_build_dir = my $tumor_build_dir = "";

		foreach my $build (@normal_builds)
		{
			my $build_status = $build->status;
			my $build_dir = $build->data_directory;
			if($build_status eq "Succeeded")
			{
				$normal_build_dir = $build_dir;
			}
		}

		foreach my $build (@tumor_builds)
		{
			my $build_status = $build->status;
			my $build_dir = $build->data_directory;
			if($build_status eq "Succeeded")
			{
				$tumor_build_dir = $build_dir;
			}
		}		
		
		my $last_build_dir = "";
		my $model_status = "New";
		my $final_build_result = "";
		my $last_build_id = 0;

		if($self->output_build_dirs || $self->output_review || $self->output_maf_file || $self->output_germline_calls || $self->uhc_filter || $self->varscan_copynumber || $self->output_loh_calls || $self->process_loh_calls)
		{
			my $num_builds = 0;		
			my $num_maf_mutations = 0;
			my $num_maf_snvs = my $num_maf_indels = 0;
	
			my $build_ids = my $build_statuses = "";
			my @builds = $model->builds;
	
			if(@builds)
			{
				$model_status = "Building";
	
				foreach my $build (@builds)
				{
					my $build_id = $build->id;
					my $build_status = $build->status;
					my $build_dir = $build->data_directory;
	
					$build_ids .= "," if($build_ids);
					$build_statuses .= "," if($build_statuses);
	
					$build_ids .= $build_id;
					$build_statuses .= $build_status;
					
					if($model_status eq "New" || $build_status eq "Succeeded" || $build_status eq "Running")
					{
						$model_status = $build_status;
						$last_build_dir = $build_dir;
					}
				}
	
				if($model->last_succeeded_build_directory)
				{
					$model_status = "Succeeded";	## Override if we have successful build dir ##				
					$succeeded_models_by_sample{$subject_name} = $model_id;
					$last_build_dir = $model->last_succeeded_build_directory;
					$last_build_dir = "/gscmnt/ams1183/info/model_data/2877873505/build113565171" if($model->id == 2877873505);
					my $tumor_model = $model->tumor_model;
					my $normal_model = $model->normal_model;
					my $tumor_sample = $tumor_model->subject_name;
					my $normal_sample = $normal_model->subject_name;

					my @temp = split(/\-/, $tumor_sample);
					my $patient_id = join("-", "TCGA", $temp[1], $temp[2]);

					if($self->output_build_dirs)
					{
						print BUILDDIRS join("\t", $tumor_sample . "-" . $normal_sample, $last_build_dir) . "\n";					
					}
	
	
					if($self->output_germline_calls)
					{
						my $dir = $self->output_germline_calls;
						if(!(-d $dir))
						{
							mkdir($dir) or die "Can't create output directory: $!\n";
						}

						## Make output dir for this model ##
						my $germline_dir = $self->output_germline_calls . "/" . $model_name;
						if(!(-d $germline_dir))
						{
							mkdir($germline_dir) or die "Can't create output directory: $!\n";
						}

						my $tumor_model = $model->tumor_model;
						my $tumor_model_dir = $tumor_model->last_succeeded_build_directory;
						my $tumor_bam = `ls $tumor_model_dir/alignments/*.bam | head -1`; chomp($tumor_bam);

						my $normal_model = $model->normal_model;
						my $normal_model_dir = $normal_model->last_succeeded_build_directory;
						my $normal_bam = `ls $normal_model_dir/alignments/*.bam | head -1`; chomp($normal_bam);

						output_germline_files($self, $last_build_dir, $germline_dir, $tumor_bam, $normal_bam, $tumor_model_dir, $normal_model_dir);
#						exit(0);
					}

					if($self->output_loh_calls)
					{
						my $dir = $self->output_loh_calls;
						if(!(-d $dir))
						{
							mkdir($dir) or die "Can't create output directory: $!\n";
						}

						## Get normal and tumor sample names ##
						my $tumor_model = $model->tumor_model;
						my $normal_model = $model->normal_model;
						my $tumor_sample = $tumor_model->subject_name;
						my $normal_sample = $normal_model->subject_name;

						## Make output dir for this model ##
						my $loh_dir = $self->output_loh_calls . "/" . $tumor_sample . "-" . $normal_sample;
						if(!(-d $loh_dir))
						{
							mkdir($loh_dir) or die "Can't create output directory: $!\n";
						}
                                                
                                                ## Determine if we already have results ##
                                                
                                                my $loh_output_file = "$loh_dir/varScan.loh.snp";
                                                if(-e $loh_output_file)
                                                {
                                                        warn "Skipping $tumor_sample-$normal_sample because file exists...\n";
                                                }
                                                else
                                                {
                                                        ## Get BAM Files ##
        
                                                        my $tumor_model_dir = $tumor_model->last_succeeded_build_directory;
                                                        my $tumor_bam = `ls $tumor_model_dir/alignments/*.bam`; chomp($tumor_bam);
        
                                                        my $normal_model_dir = $normal_model->last_succeeded_build_directory;
                                                        my $normal_bam = `ls $normal_model_dir/alignments/*.bam`; chomp($normal_bam);
        
                                                        ## Output LOH Files ##
                                                        
                                                        output_loh_files($self, $last_build_dir, $loh_dir, $normal_bam, $tumor_bam);                                                        
                                                }


#						$debug_counter++;
#						exit(0) if($debug_counter >= 10);
					}
					elsif($self->process_loh_calls)
					{
						## Get normal and tumor sample names ##
						my $tumor_model = $model->tumor_model;
						my $normal_model = $model->normal_model;
						my $tumor_sample = $tumor_model->subject_name;
						my $normal_sample = $normal_model->subject_name;

						my $loh_dir = $self->process_loh_calls . "/" . $tumor_sample . "-" . $normal_sample;
						
						if(-d $loh_dir)
						{
                                                        my $merged_output_file = "$loh_dir/varScan.merged.snp";
                                                        if(-e $merged_output_file)
                                                        {
                                                                warn "Skipping $tumor_sample-$normal_sample because merged file exists...\n";       
                                                        }
                                                        else
                                                        {
        							process_loh($loh_dir);
                						sleep(1);
                                                        }
                                                }

					}
					
					my %build_results = get_build_results($last_build_dir);
					$final_build_result = $build_results{'tier1_snvs'} . " Tier1 SNVs, " . $build_results{'tier1_indels'} . " Tier1 Indels, ";
	
					if($self->output_review)
					{
						my $tumor_model = $model->tumor_model;
						my $normal_model = $model->normal_model;
						my $tumor_model_dir = $tumor_model->last_succeeded_build_directory;
						my $normal_model_dir = $normal_model->last_succeeded_build_directory;
						my $tumor_bam = `ls $tumor_model_dir/alignments/*.bam`; chomp($tumor_bam);
						my $normal_bam = `ls $normal_model_dir/alignments/*.bam`; chomp($normal_bam);
						
						my $tier1_snvs = $last_build_dir . "/merged.somatic.snp.filter.novel.tier1";

						## If UHC-filter flag turned on and file exists, use it instead ##
						if($self->uhc_filter && -e "$tier1_snvs.uhc-filter")
						{
							my $uhc_count = `cat $tier1_snvs.uhc-filter | wc -l`; chomp($uhc_count);
							$stats{'uhc_sites_skipped'} += $uhc_count;
							$tier1_snvs = $last_build_dir . "/merged.somatic.snp.filter.novel.tier1.uhc-filter.removed";
						}
						elsif($self->uhc_filter)
						{
							die "Warning: $patient_id has no UHC-filter file in $last_build_dir\n";
						}

						my $output_tier1_snvs = $self->output_review . "/" . $subject_name . ".$model_id.SNVs.tsv";
						output_snvs_for_review($self, $model_id, $tier1_snvs, $output_tier1_snvs, $subject_name, $normal_bam, $tumor_bam);
	
						my $tier1_gatk = $last_build_dir . "/gatk.output.indel.formatted.Somatic.tier1";
						my $tier1_indels = $last_build_dir . "/merged.somatic.indel.filter.tier1";
						## Grab HC Indels instead if possible ##
						$tier1_indels .= ".hc" if(-e "$tier1_indels.hc"); 
						
						my $output_tier1_indels = $self->output_review . "/" . $subject_name . ".$model_id.Indels.tsv";
						output_indels_for_review($model_id, $tier1_indels, $tier1_gatk, $output_tier1_indels, $subject_name, $normal_bam, $tumor_bam);
					}
					elsif($self->uhc_filter)
					{
						my $tumor_model = $model->tumor_model;
						my $normal_model = $model->normal_model;
						my $tumor_model_dir = $tumor_model->last_succeeded_build_directory;
						my $normal_model_dir = $normal_model->last_succeeded_build_directory;
						my $tumor_bam = `ls $tumor_model_dir/alignments/*.bam`; chomp($tumor_bam);
						my $normal_bam = `ls $normal_model_dir/alignments/*.bam`; chomp($normal_bam);
						
						my $tier1_snvs = $last_build_dir . "/merged.somatic.snp.filter.novel.tier1";
						my $output_file = $tier1_snvs . ".uhc-filter";
						my $filtered_file = $tier1_snvs . ".uhc-filter.removed";
						
						## Save UHC SNVs or else generate them ##
						if(-e $output_file)
						{
							save_uhc_calls($patient_id, $output_file);
						}
						else
						{
							if($self->reference)
							{
								system("bsub -q short -R\"select[model!=Opteron250 && mem>6000] rusage[mem=6000]\" -M 6000000 gmt somatic ultra-high-confidence --min-tumor-var-freq 0.10 --tumor-bam $tumor_bam --normal-bam $normal_bam --variant-file $tier1_snvs --output-file $output_file --filtered-file $filtered_file --reference " . $self->reference);
							}
							else
							{
								system("bsub -q short -R\"select[model!=Opteron250 && mem>6000] rusage[mem=6000]\" -M 6000000 gmt somatic ultra-high-confidence --min-tumor-var-freq 0.10 --tumor-bam $tumor_bam --normal-bam $normal_bam --variant-file $tier1_snvs --output-file $output_file --filtered-file $filtered_file");															
							}

						}
					}

					if($self->uhc_indels)
					{
						my $tier1_indels = $last_build_dir . "/merged.somatic.indel.filter.tier1";
						if(-e "$tier1_indels.hc")
						{
							save_uhc_calls($patient_id, "$tier1_indels.hc");
						}
					}

					if($self->varscan_copynumber)
					{
						my $tumor_model = $model->tumor_model;
						my $normal_model = $model->normal_model;
						my $tumor_sample = $tumor_model->subject_name;
						my $normal_sample = $normal_model->subject_name;
						my $tumor_model_dir = $tumor_model->last_succeeded_build_directory;
						my $normal_model_dir = $normal_model->last_succeeded_build_directory;
						my $tumor_bam = `ls $tumor_model_dir/alignments/*.bam`; chomp($tumor_bam);
						my $normal_bam = `ls $normal_model_dir/alignments/*.bam`; chomp($normal_bam);
					
						my $output_dir = $self->varscan_copynumber . "/" . $tumor_sample . "-" . $normal_sample;
						mkdir($output_dir) if(!(-d $output_dir));
						my $cmd = "gmt varscan copy-number --output $output_dir/varScan.output --normal-bam $normal_bam --tumor-bam $tumor_bam";
						$cmd .= " --reference " . $self->reference if($self->reference);
						if(!(-e "$output_dir/varScan.output.copynumber"))
						{
							system("bsub -q long -R\"select[model!=Opteron250 && mem>4000]\" $cmd");							
						}
						else
						{
							print "Skipping $tumor_sample-$normal_sample...\n";
						}


					}
					
					if($self->output_maf_file)
					{
						my $sample_maf_file = $last_build_dir . "/tcga-maf.tsv";
						
						## If no MAF file exists, let's build one ##
						if(!(-e $sample_maf_file))
						{
							print "MAF File not found, so generating one...\n";
							my $cmd_obj = Genome::Model::Tools::Capture::BuildMafFile->create(
								data_dir => , $last_build_dir,
								normal_sample => $normal_subject_name,
								tumor_sample => $tumor_subject_name,					    
								output_file => $sample_maf_file,
							);
	
							$cmd_obj->execute;	
						}
						
						if(-e $sample_maf_file)
						{
							my $sample_maf_results = parse_maf_file($subject_name, $sample_maf_file);
							my @sample_results = split(/\n/, $sample_maf_results);
							$num_maf_mutations = @sample_results;
							
							foreach my $maf_line (@sample_results)
							{
								my @temp = split(/\t/, $maf_line);
								my $variant_type = $temp[9];
								$num_maf_snvs++ if($variant_type eq "SNP" || $variant_type eq "SNV");
								$num_maf_indels++ if($variant_type eq "INS" || $variant_type eq "DEL");
							}
							
							if($maf_header && !$maf_header_printed)
							{
								print MAF "$maf_header\n";
								$maf_header_printed = 1;
							}
							print MAF "$sample_maf_results\n";
						}
	
					}
	
					
				}
	
			}
	
			## Figure out how many sites have already been reviewed ##
			my $num_sites_passed_review = $patient_sites_passed_review{$patient_id};
			$num_sites_passed_review = 0 if(!$num_sites_passed_review);
		}

		print join("\t", $model_id, $model_name, $normal_model_id, $tumor_model_id, $final_build_result) . "\n";

#		print join("\t", $model_id, $subject_name, $model_status, $build_ids, $build_statuses, $final_build_result, $num_maf_mutations . " mutations added to MAF", $num_maf_snvs . " SNVs", $num_maf_indels . " Indels") . "\n";
#		print join("\t", $model_name, $model_id, $build_ids, $build_statuses, $final_build_result, $num_maf_mutations, $num_maf_snvs, $num_maf_indels) . "\n";
#		print join("\t", $model_name, $model_id, $final_build_result, "$num_sites_passed_review pass review", "$num_maf_mutations mutations ($num_maf_snvs SNVs, $num_maf_indels indels) added to MAF") . "\n";
#		print "$normal_build_dir\n$tumor_build_dir\n";

	}	
	
	print $stats{'models_in_group'} . " models in group\n" if($stats{'models_in_group'});
	print $stats{'models_running'} . " models running\n" if($stats{'models_running'});
	print $stats{'models_finished'} . " models finished\n" if($stats{'models_finished'});

	if($self->output_review)
	{
		print $stats{'review_snvs_possible'} . " Tier 1 SNVs could be reviewed\n";
		print $stats{'uhc_sites_skipped'} . " UHC SNVs were excluded from review\n";
		print $stats{'review_snvs_already'} . " were already reviewed\n";
		print $stats{'review_snvs_already_wildtype'} . " were wild-type in another sample\n";
		print $stats{'review_snvs_already_germline'} . " were germline in at least 3 other samples\n";
		print $stats{'review_snvs_filtered'} . " were filtered as probable germline\n";
		print $stats{'review_snvs_uhc_lowcov'} . " were removed because low-coverage makes them ambiguous\n";
		print $stats{'review_snvs_uhc_normalfreq'} . " were removed due to variant presence in normal\n";
		print $stats{'review_snvs_uhc_freqdiff'} . " were removed due to low freq diff between tumor and normal\n";
		
		print $stats{'review_snvs_included'} . " were included for review\n";
		print $stats{'review_snvs_already_included'} . " were duplicates and not counted twice\n";

		print $stats{'review_indels_possible'} . " Tier 1 Indels could be reviewed\n";
		print $stats{'review_indels_already'} . " were already reviewed\n";
		print $stats{'review_indels_already_wildtype'} . " were wild-type in another sample\n";
		print $stats{'review_indels_already_germline'} . " were germline in another sample\n";		
		print $stats{'review_indels_filtered'} . " were filtered as probable germline\n";
		print $stats{'review_indels_included'} . " were included for review\n";
	}


	close(MAF) if($self->output_maf_file);
	
	return 1;
}



################################################################################################
# Get Build Results - Summarize the progress/results of a given build
#
################################################################################################

sub save_uhc_calls
{
	my $patient_id = shift(@_);
	my $FileName = shift(@_);

	## Parse the Tier 1 SNVs file ##

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($chrom, $chr_start, $chr_stop) = split(/\t/, $line);
		my $variant_key = join("\t", $patient_id, $chrom, $chr_start, $chr_stop);
		$passed_sites{$variant_key} = 1;
	}
	
	close($input);
}


################################################################################################
# Get Build Results - Summarize the progress/results of a given build
#
################################################################################################

sub parse_maf_file
{
	my $sample_name = shift(@_);
	my $FileName = shift(@_);

	my $sample_maf = "";

	## Column index for fields in MAF file ##
	
	my %column_index = ();
	my @columns = ();

	## Parse the Tier 1 SNVs file ##

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		my @lineContents = split(/\t/, $line);
	
		if($lineCounter <= 2 && $line =~ "Chrom")
		{
			$maf_header = $line;
			
			## Parse the MAF header line to determine field locations ##	
			my $numContents = @lineContents;
			
			for(my $colCounter = 0; $colCounter < $numContents; $colCounter++)
			{
				if($lineContents[$colCounter])
				{
					$column_index{$lineContents[$colCounter]} = $colCounter;
				}
			}
			
			foreach my $column (keys %column_index)
			{
				## Print out the columns as parsed ##
				#print "$column_index{$column}\t$column\n";
				$columns[$column_index{$column}] = $column;	## Save the column order ##
			}
		}
		elsif($lineCounter < 2)
		{

		}
		elsif($lineCounter > 2 && !@columns)
		{
			die "No Header in MAF file!\n";
		}
		elsif($lineCounter > 2 && @columns)
		{
			## Build a record for this line, assigning all values to respective fields ##
			
			my %record = ();

			foreach my $column_name (keys %column_index)
			{
				my $column_number = $column_index{$column_name};
				$record{$column_name} = $lineContents[$column_number];
			}
			
			my $chrom = $record{'Chromosome'};
			my $chr_start = $record{'Start_position'};
			my $chr_stop = $record{'End_position'};
			my $ref_allele = $record{'Reference_Allele'};
			my $var_allele = $record{'Tumor_Seq_Allele2'};
			$var_allele = $record{'Tumor_Seq_Allele1'} if($var_allele eq $ref_allele);
			my $var_type = $record{'Variant_Type'};
			
			my @temp = split(/\-/, $sample_name);
			my $patient_id = join("-", "TCGA", $temp[1], $temp[2]);
			my $variant_key = join("\t", $patient_id, $chrom, $chr_start, $chr_stop); #, $ref_allele, $var_allele);
			
			## Include variant if it had a review-passed call, or if no reviews were loaded ##
			if($passed_sites{$variant_key} || !(%passed_sites))
			{
				$sample_maf .= "\n" if($sample_maf);
				$sample_maf .= $line;
			}
		}
	}
	
	close($input);
	
	return($sample_maf);
}


################################################################################################
# Get Build Results - Summarize the progress/results of a given build
#
################################################################################################

sub get_build_results
{
	my $build_dir = shift(@_);
	my %results = ();
	
	my $tier1_snvs = $build_dir . "/merged.somatic.snp.filter.novel.tier1";
	my $tier1_gatk = $build_dir . "/gatk.output.indel.formatted.Somatic.tier1";
	my $tier1_indels = $build_dir . "/merged.somatic.indel.filter.tier1";
	
	## Check for Tier 1 SNVs ##
	
	if(-e $tier1_snvs)
	{
		## Get count ##
		my $count = `cat $tier1_snvs | wc -l`;
		chomp($count);
		$results{'tier1_snvs'} = $count;
	}
	
	if(-e $tier1_indels && -e $tier1_gatk)
	{
		my $count = `cat $tier1_indels $tier1_gatk | cut --fields=1-3 | sort -u | wc -l`;
		chomp($count);

		$results{'tier1_indels'} = $count;
	}

	if(-e $tier1_indels)
	{
		my $count = `cat $tier1_indels | cut --fields=1-3 | sort -u | wc -l`;
		chomp($count);

		$results{'tier1_indels'} = $count;
	}
	
	return(%results);
}





################################################################################################
# Get Build Results - Summarize the progress/results of a given build
#
################################################################################################

sub get_build_progress
{
	my $build_dir = shift(@_);
	my %results = ();
	
	my $tier1_snvs = $build_dir . "/merged.somatic.snp.filter.novel.tier1";
	my $tier1_gatk = $build_dir . "/gatk.output.indel.formatted.Somatic.tier1";
	my $tier1_indels = $build_dir . "/merged.somatic.indel.filter.tier1";
	
	## Check for Tier 1 SNVs ##
	
	if(-e $tier1_snvs)
	{
		$results{'tier1_snvs_done'} = 1;
		## Get count ##
		my $count = `cat $tier1_snvs | wc -l`;
		chomp($count);
		$results{'num_tier1_snvs'} = $count;
	}
	
	return(%results);
}


################################################################################################
# Get Build Results - Summarize the progress/results of a given build
#
################################################################################################

sub load_review_database
{
	my $FileName = shift(@_);
	my $num_passed_sites = 0;
	## Parse the Tier 1 SNVs file ##

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		my ($model_id, $build_id, $sample_name, $chrom, $chr_start, $chr_stop, $ref, $var, $code) = split(/\t/, $line);
		
		## Determine if this is a BED-formatted file ##
		
		if($chr_start ne "chr_start" && $chr_stop - $chr_start == 1)
		{
			## Especially, for SNVs ##
			if(length($ref) == 1 && length($var) == 1 && $ref ne "-" && $ref ne "0" && $var ne "-" && $var ne "0")
			{
				$chr_start++;
			}
		}
		
		## Abbreviate sample name ##
		
		my @temp = split(/\-/, $sample_name);
		my $patient_id = join("-", "TCGA", $temp[1], $temp[2]);
		
		my $key = join("\t", $patient_id, $chrom, $chr_start, $chr_stop);
#		my $key = join("\t", $model_id, $chrom, $chr_start, $chr_stop);
		
		if($code ne "A")
		{
			$already_reviewed{$key} = $code;			
		}

		if($code eq "O" || $code eq "D" || $code eq "LQ")
		{
			my $key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
			$wildtype_sites{$key}++;
		}
		elsif($code eq "G")
		{
			my $key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
			$germline_sites{$key}++;			
		}
		elsif($code eq "S" || $code eq "V")
		{
			my @temp = split(/\-/, $sample_name);
			my $patient_id = join("-", "TCGA", $temp[1], $temp[2]);
			if(!$temp[1] || !$temp[2])
			{
				die "Error parsing line from $FileName: $line had incomplete sample name $sample_name\n";
			}
			my $key = join("\t", $patient_id, $chrom, $chr_start, $chr_stop);
			$passed_sites{$key} = 1;
			$num_passed_sites++;
			
			$patient_sites_passed_review{$patient_id}++;
			
		}
	}
	
	close($input);
	
	print "$num_passed_sites sites passed review\n";
}



################################################################################################
# Get Build Results - Summarize the progress/results of a given build
#
################################################################################################

sub output_snvs_for_review
{
	my ($self, $model_id, $variant_file, $output_file, $subject_name, $normal_bam, $tumor_bam) = @_;
	
	## Check for Tier 1 SNVs ##
	
	if(-e $variant_file)
	{
		## Open the output file ##
		
		open(OUTFILE, ">$output_file") or die "Can't open output file: $!\n";
		print OUTFILE join("\t", "TUMOR", $tumor_bam) . "\n";
		print OUTFILE join("\t", "NORMAL", $normal_bam) . "\n";
		print OUTFILE "chrom\tchr_start\tchr_stop\tref\tvar\tcode\tnote\n";
		
		## Parse the Tier 1 SNVs file ##
	
		my $input = new FileHandle ($variant_file);
		my $lineCounter = 0;
	
		while (<$input>)
		{
			chomp;
			my $line = $_;
			$lineCounter++;

			my ($chrom, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $line);
			my @lineContents = split(/\t/, $line);

			my $include_flag = 0;

			## Get Patient ID ##

			my @temp = split(/\-/, $subject_name);
			my $patient_id = join("-", "TCGA", $temp[1], $temp[2]);
			
			my $key = join("\t", $patient_id, $chrom, $chr_start, $chr_stop);		
#			my $key = join("\t", $model_id, $chrom, $chr_start, $chr_stop);
			my $variant_key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
			my $sample_variant_key = join("\t", $subject_name, $variant_key);
			
			$stats{'review_snvs_possible'}++;
			
			if($already_reviewed{$key})
			{
				$include_flag = 0;
				$stats{'review_snvs_already'}++;
			}
			elsif($wildtype_sites{$variant_key})
			{
				$include_flag = 0;
				$stats{'review_snvs_already_wildtype'}++;
			}
			elsif($germline_sites{$variant_key} && $germline_sites{$variant_key} >= 3)
			{
				$include_flag = 0;
				$stats{'review_snvs_already_germline'}++;
			}
			elsif($included_variants{$sample_variant_key})
			{
				$include_flag = 0;
				$stats{'review_snvs_already_included'}++;				
			}
			else
			{
				## Sniper SNV/INS/DEL File ##
				
				if($lineContents[5] && ($lineContents[5] eq "SNP" || $lineContents[5] eq "INS" || $lineContents[5] eq "DEL"))
				{
					$include_flag = 1;
				}
				
				## VarScan File ##
				
				elsif($line =~ 'Somatic')
				{
					my $normal_freq = $lineContents[7];
					my $tumor_freq = $lineContents[11];
					$normal_freq =~ s/\%//g;
					$tumor_freq =~ s/\%//g;
					
					if($tumor_freq < 30 && $normal_freq >= 2)
					{
						## Exclude a possible Germline Event ##
						$stats{'review_snvs_filtered'}++;
					}
					else
					{
						$include_flag = 1;			
					}
				}
				
				## GATK Indel File ##
				
				elsif($line =~ 'OBS\_COUNTS')
				{
					$include_flag = 1;				
				}				
			}
			
			## If the UHC filter was applied, do further processing ##
			
			if($self->uhc_filter && $include_flag)
			{
				my @lineContents = split(/\t/, $line);
				my $numContents = @lineContents;
				
				for(my $colCounter = 5; $colCounter < $numContents; $colCounter++)
				{
					## IF failed due to coverage, handle it ##
					if($lineContents[$colCounter] && $lineContents[$colCounter] eq "Failed:Coverage")
					{
						my $normal_cov = $lineContents[$colCounter + 1];
						my $tumor_cov = $lineContents[$colCounter + 2];
						
						if($normal_cov < 4 || $tumor_cov < 6)
						{
							$include_flag = 0;
							$stats{'review_snvs_uhc_lowcov'}++;
						}
					}
					elsif($lineContents[$colCounter] && $lineContents[$colCounter] eq "Failed:VarFreq")
					{
						my $normal_cov = $lineContents[$colCounter + 1];
						my $tumor_cov = $lineContents[$colCounter + 2];
						my $normal_freq = $lineContents[$colCounter + 3];
						my $tumor_freq = $lineContents[$colCounter + 4];
						$normal_freq =~ s/\%//;
						$tumor_freq =~ s/\%//;
						my $freq_diff = $tumor_freq - $normal_freq;
						if($normal_freq > 5 && $normal_cov >= 6)
						{
							$include_flag = 0;
							$stats{'review_snvs_uhc_normalfreq'}++;
						}
						elsif($freq_diff < 10)
						{
							$include_flag = 0;
							$stats{'review_snvs_uhc_freqdiff'}++;
						}
					}
				}
			}
			
			

			if($include_flag)
			{
				print OUTFILE join("\t", $chrom, $chr_start, $chr_stop, $ref, $var) . "\n";
				$stats{'review_snvs_included'}++;
				my $key = join("\t", $subject_name, $chrom, $chr_start, $chr_stop, $ref, $var);
				$included_variants{$sample_variant_key} = 1;
			}


		}
		
		close($input);
		
		
		close(OUTFILE);
		
	}

}



################################################################################################
# Get Build Results - Summarize the progress/results of a given build
#
################################################################################################

sub output_indels_for_review
{
	my ($model_id, $variant_file1, $variant_file2, $output_file, $subject_name, $normal_bam, $tumor_bam) = @_;
	
	my %indels = ();
	
	## Check for Tier 1 SNVs ##
	
	if(-e $variant_file1)
	{
		## Parse the Tier 1 SNVs file ##
	
		my $input = new FileHandle ($variant_file1);
		my $lineCounter = 0;
	
		while (<$input>)
		{
			chomp;
			my $line = $_;
			$lineCounter++;

			my ($chrom, $chr_start) = split(/\t/, $line);
		
			$indels{"$chrom\t$chr_start"} .= "\n" if($indels{"$chrom\t$chr_start"});
			$indels{"$chrom\t$chr_start"} .= $line;
		}
		
		close($input);		
	}


	if($variant_file2 && -e $variant_file2)
	{
		## Parse the Tier 1 SNVs file ##
	
		my $input = new FileHandle ($variant_file2);
		my $lineCounter = 0;
	
		while (<$input>)
		{
			chomp;
			my $line = $_;
			$lineCounter++;

			my ($chrom, $chr_start) = split(/\t/, $line);
		
			$indels{"$chrom\t$chr_start"} .= "\n" if($indels{"$chrom\t$chr_start"});
			$indels{"$chrom\t$chr_start"} .= $line;
		}
		
		close($input);		
	}


	## Open the output file ##
	
	open(OUTFILE, ">$output_file") or die "Can't open output file: $!\n";
	print OUTFILE join("\t", "TUMOR", $tumor_bam) . "\n";
	print OUTFILE join("\t", "NORMAL", $normal_bam) . "\n";
	print OUTFILE "chrom\tchr_start\tchr_stop\tref\tvar\tcode\tnote\n";

	foreach my $key (sort byChrPos keys %indels)
	{
		my $include_flag = 0;
		
		my ($chrom, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $indels{$key});
		
		my $review_key = join("\t", $model_id, $chrom, $chr_start, $chr_stop);
		my $variant_key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
		
		my @indel_lines = split(/\n/, $indels{$key});
		my $num_indel_lines = @indel_lines;
		
		$stats{'review_indels_possible'}++;
		
		if($already_reviewed{$review_key})
		{
			## Skip already reviewed ##
			$stats{'review_indels_already'}++;
			$include_flag = 0;
		}
		elsif($wildtype_sites{$variant_key})
		{
			$include_flag = 0;
			$stats{'review_indels_already_wildtype'}++;
		}
		elsif($germline_sites{$variant_key}) # && $germline_sites{$variant_key} >= 3)
		{
			$include_flag = 0;
			$stats{'review_indels_already_germline'}++;
		}		
		elsif($indels{$key} =~ 'OBS\_COUNTS')
		{
			## GATK indel ##
			$include_flag = 1;
		}
		else
		{
			my $line = $indels{$key};
			my @lineContents = split(/\t/, $line);			

			## Sniper SNV/INS/DEL File ##
			
			if($lineContents[5] && ($lineContents[5] eq "SNP" || $lineContents[5] eq "INS" || $lineContents[5] eq "DEL"))
			{
				$include_flag = 1;
			}
			
			## VarScan File ##
			
			elsif($line =~ 'Somatic')
			{
				my $normal_freq = $lineContents[7];
				my $tumor_freq = $lineContents[11];
				$normal_freq =~ s/\%//g;
				$tumor_freq =~ s/\%//g;
				
				if($tumor_freq < 20 && $normal_freq >= 2)
				{
					$stats{'review_indels_filtered'}++;
				}
				else
				{
					$include_flag = 1;
				}
			}
		}
		
		if($include_flag)
		{
			print OUTFILE join("\t", $chrom, $chr_start, $chr_stop, $ref, $var) . "\n";
			$stats{'review_indels_included'}++;
		}
	}
	
	close(OUTFILE);


}



################################################################################################
# Output Germline Files - output the files for germline analysis 
# 1.) Get BED files from SAMtools -> Subtract Somatic -> FPfilter -> Combine
# 2.) Get Germline/LOH calls from VarScan -> Combine -> FPfilter
################################################################################################

sub output_germline_files
{
	my ($self, $build_dir, $germline_dir, $tumor_bam, $normal_bam, $tumor_model_dir, $normal_model_dir) = @_;

	my $cmd = "";	
#	print "Somatic\t$build_dir\nTumor\t$tumor_model_dir\nNormal\t$normal_model_dir\n";

	## Open output script ##
	
	open(SCRIPT, ">$germline_dir/germline.sh");
	print SCRIPT "#!/gsc/bin/sh\n";

	## Get File of Somatic SNVs ##
	
	my $somatic_snvs = "$build_dir/merged.somatic.snp.filter.novel";

	## Get SAMtools variant files with Annotation ##
	my $normal_annotation_file = "$normal_model_dir/variants/filtered.variants.post_annotation";
	$normal_annotation_file = "$normal_model_dir/snp_related_metrics/filtered.variants.post_annotation" if(!(-e $normal_annotation_file));

	my $tumor_annotation_file = "$tumor_model_dir/variants/filtered.variants.post_annotation";
	$tumor_annotation_file = "$tumor_model_dir/snp_related_metrics/filtered.variants.post_annotation" if(!(-e $tumor_annotation_file));

	## Get SAMtools variant files with Annotation ##
	my $normal_variant_file = "$normal_model_dir/variants/snps_all_sequences.filtered";
	$normal_variant_file = "$normal_model_dir/snp_related_metrics/snps_all_sequences.filtered" if(!(-e $normal_variant_file));
	$normal_variant_file = "$normal_model_dir/variants/snv/samtools-r599-/snp-filter-v1-/snvs.hq" if(!(-e $normal_variant_file));

	my $tumor_variant_file = "$tumor_model_dir/variants/snps_all_sequences.filtered";
	$tumor_variant_file = "$tumor_model_dir/snp_related_metrics/snps_all_sequences.filtered" if(!(-e $tumor_variant_file));
	$tumor_variant_file = "$tumor_model_dir/variants/snv/samtools-r599-/snp-filter-v1-/snvs.hq" if(!(-e $tumor_variant_file));
	
	## Format Normal and Tumor Variant Files ##
	print SCRIPT "echo Isolating normal/tumor SNP calls...\n";
	$cmd = "gmt capture format-snvs --variant $normal_variant_file --output-file $germline_dir/samtools.normal.snp --preserve-call 1 --append-line 0";
	print SCRIPT "$cmd\n";
	$cmd = "gmt capture format-snvs --variant $tumor_variant_file --output-file $germline_dir/samtools.tumor.snp --preserve-call 1 --append-line 0";
	print SCRIPT "$cmd\n";
	
	## Get VarScan Files ##

	print SCRIPT "echo Getting VarScan files...\n";
	my $varscan_germline_snps = `ls $build_dir/varScan.output.snp.formatted.Germline.hc 2>/dev/null`;
	chomp($varscan_germline_snps);

	my $varscan_loh_snps = `ls $build_dir/varScan.output.snp.formatted.LOH.hc 2>/dev/null`;
	chomp($varscan_loh_snps);

	$cmd = "cut -f 1-4,9 $varscan_germline_snps >$germline_dir/varScan.germline.snp";
	print SCRIPT "$cmd\n";
	$cmd = "cut -f 1-4,9 $varscan_loh_snps >$germline_dir/varScan.loh.snp";
	print SCRIPT "$cmd\n";

	## Merge files ##
	
	print SCRIPT "echo Merging files...\n";
	$cmd = "cat $germline_dir/samtools.normal.snp $germline_dir/varScan.germline.snp $germline_dir/varScan.loh.snp | cut -f 1-5 | sort -u >$germline_dir/merged.normal.snp";
	print SCRIPT "$cmd\n";
	$cmd = "cat $germline_dir/samtools.tumor.snp $germline_dir/varScan.germline.snp | cut -f 1-5 | sort -u >$germline_dir/merged.tumor.snp";
	print SCRIPT "$cmd\n";

	## Sort merged files ##

	print SCRIPT "echo Sorting files...\n";	
	$cmd = "gmt capture sort-by-chr-pos --input $germline_dir/merged.normal.snp --output $germline_dir/merged.normal.snp";
	print SCRIPT "$cmd\n";
	$cmd = "gmt capture sort-by-chr-pos --input $germline_dir/merged.tumor.snp --output $germline_dir/merged.tumor.snp";
	print SCRIPT "$cmd\n";

	## Convert to bed ##
	
	print SCRIPT "echo Converting to BED...\n";
	$cmd = "gmt capture convert-to-bed --variant $germline_dir/merged.normal.snp --output $germline_dir/merged.normal.snp.bed";
	print SCRIPT "$cmd\n";
	$cmd = "gmt capture convert-to-bed --variant $germline_dir/merged.tumor.snp --output $germline_dir/merged.tumor.snp.bed";
	print SCRIPT "$cmd\n";

	print SCRIPT "echo Fast-tiering...\n";
	$cmd = "gmt fast-tier fast-tier --variant-bed-file $germline_dir/merged.normal.snp.bed --tier1-output $germline_dir/merged.normal.snp.tier1.bed --tier-file-location " . $self->tier_file_location;
	print SCRIPT "$cmd\n";;

	$cmd = "gmt fast-tier fast-tier --variant-bed-file $germline_dir/merged.tumor.snp.bed --tier1-output $germline_dir/merged.tumor.snp.tier1.bed --tier-file-location " . $self->tier_file_location;
	print SCRIPT "$cmd\n";

	print SCRIPT "echo Converting from BED...\n";
	$cmd = "gmt capture convert-from-bed --variant-file $germline_dir/merged.normal.snp.tier1.bed --output-file $germline_dir/merged.normal.snp.tier1";
	print SCRIPT "$cmd\n";
	$cmd = "gmt capture convert-from-bed --variant-file $germline_dir/merged.tumor.snp.tier1.bed --output-file $germline_dir/merged.tumor.snp.tier1";
	print SCRIPT "$cmd\n";

	print SCRIPT "echo Running FP Filter...\n";
	## Apply FP-filter to Germline using Tumor BAM ##
	#bsub -q long -R\"select[type==LINUX64 && mem>8000] rusage[mem=8000]\" -M 8000000 
	$cmd = "gmt somatic filter-false-positives --reference " . $self->reference . " --max-mm-qualsum-diff 100 ";
	$cmd .= "--variant-file $germline_dir/merged.tumor.snp.tier1 --bam-file $tumor_bam --output-file $germline_dir/merged.tumor.snp.tier1.fpfilter --filtered-file $germline_dir/merged.tumor.snp.tier1.fpfilter.removed";
	print SCRIPT "$cmd\n";

	$cmd = "gmt somatic filter-false-positives --reference " . $self->reference . " --max-mm-qualsum-diff 100 ";
	$cmd .= "--variant-file $germline_dir/merged.normal.snp.tier1 --bam-file $normal_bam --output-file $germline_dir/merged.normal.snp.tier1.fpfilter --filtered-file $germline_dir/merged.normal.snp.tier1.fpfilter.removed";
	print SCRIPT "$cmd\n";

	$cmd = "gmt capture combine-snv-files --variant-file1 $germline_dir/merged.normal.snp.tier1.fpfilter --variant-file2 $germline_dir/merged.tumor.snp.tier1.fpfilter --output-file $germline_dir/merged.germline.snp";
	print SCRIPT "$cmd\n";

	$cmd = "gmt annotate transcript-variants --annotation-filter top --reference-transcripts " . $self->reference_transcripts . " ";
	$cmd .= "--variant-file $germline_dir/merged.germline.snp --output-file $germline_dir/merged.germline.snp.transcript-annotation";
	print SCRIPT "$cmd\n";

	## Merge SNVs with Annotation ##
	print SCRIPT "echo Merging with Annotation...\n";	
	$cmd = "gmt analysis somatic-pipeline merge-snvs-with-annotation --variant $germline_dir/merged.germline.snp --output $germline_dir/merged.germline.snp.annotated --annotation $germline_dir/merged.germline.snp.transcript-annotation";
	print SCRIPT "$cmd\n";

	close(SCRIPT);
	
	system("bsub -q long -R\"select[type==LINUX64 && mem>8000] rusage[mem=8000]\" -M 8000000 -oo $germline_dir/germline.sh.log \"sh $germline_dir/germline.sh\"");


	## Process Germline Indels ##

	open(SCRIPT, ">$germline_dir/germline-indel.sh");
	print SCRIPT "#!/gsc/bin/sh\n";

	## Get GATK File ##
	
	my $gatk_indel_file = `ls $build_dir/gatk.output.indel.formatted`;
	chomp($gatk_indel_file);

	my $gatk_output_file = "$germline_dir/gatk.germline.indel";
	print SCRIPT "echo Getting germline GATK calls...\n";
	if($gatk_indel_file && !(-e $gatk_output_file))
	{	
		if($self->germline_roi_file)
		{
			$cmd = "grep GERMLINE $gatk_indel_file >$gatk_output_file.temp";
			print SCRIPT "$cmd\n";
			print SCRIPT "echo Limiting to ROI...\n";
			$cmd = "java -jar $ENV{GENOME_SW_LEGACY_JAVA}/VarScan/VarScan.v2.2.6.jar limit $gatk_output_file.temp --regions-file " . $self->germline_roi_file . " --output-file $gatk_output_file";
			print SCRIPT "$cmd\n";
			$cmd = "rm -rf $gatk_output_file.temp";
			print SCRIPT "$cmd\n";
		}
		else
		{
			## Copy ALL Indels ##
			$cmd = "grep GERMLINE $gatk_indel_file >$gatk_output_file";
			print SCRIPT "$cmd\n";
		}

	}

	my $varscan_germline_indels = `ls $build_dir/varScan.output.indel.formatted.Germline.hc 2>/dev/null`;
	chomp($varscan_germline_indels);

	my $varscan_loh_indels = `ls $build_dir/varScan.output.indel.formatted.LOH.hc 2>/dev/null`;
	chomp($varscan_loh_indels);
	print SCRIPT "echo Getting VarScan indel calls...\n";
	if(-e $varscan_germline_indels && -e $varscan_loh_indels)
	{
		my $output_indel_varscan = $germline_dir . "/varScan.germline.indel";
		$cmd = "cat $varscan_germline_indels $varscan_loh_indels >$output_indel_varscan";
		print SCRIPT "$cmd\n";
		$cmd = "gmt capture sort-by-chr-pos --input $output_indel_varscan --output $output_indel_varscan";
		print SCRIPT "$cmd\n";

		if($self->germline_roi_file)
		{
			print SCRIPT "echo Limiting to ROI...\n";
			$cmd = "java -jar $ENV{GENOME_SW_LEGACY_JAVA}/VarScan/VarScan.jar limit $output_indel_varscan --regions-file " . $self->germline_roi_file . " --output-file $output_indel_varscan.roi";
			print SCRIPT "$cmd\n";
			$output_indel_varscan .= ".roi";
		}
		
		print SCRIPT "echo Filtering false indels...\n";
		## Apply FP-filter using Normal BAM ##
		$cmd = "gmt somatic filter-false-indels ";
		$cmd .= "--reference " . $self->reference . " ";
		$cmd .= "--variant-file $output_indel_varscan --bam-file $normal_bam --output-file $output_indel_varscan.fpfilter --filtered-file $output_indel_varscan.fpfilter.removed ";
		$cmd .= "--max-mm-qualsum-diff 100 ";
		print SCRIPT "$cmd\n";

		print SCRIPT "echo Combining Indel Files...\n";
		$cmd = "gmt capture combine-indel-files --variant-file1 $output_indel_varscan.fpfilter --variant-file2 $gatk_output_file --output-file $germline_dir/merged.germline.indel";
		print SCRIPT "$cmd\n";

		## Also run transcript-annotation ##
		print SCRIPT "echo Running Transcript Annotation...\n";
		$cmd = "gmt annotate transcript-variants --annotation-filter top --reference-transcripts " . $self->reference_transcripts . " ";
		$cmd .= "--variant-file $germline_dir/merged.germline.indel --output-file $germline_dir/merged.germline.indel.transcript-annotation";
		print SCRIPT "$cmd\n";

		## Merge indels with annotation ##
		print SCRIPT "echo Merging Indels with Annotation...\n";
		$cmd = "gmt analysis somatic-pipeline merge-indels-with-annotation --variant $germline_dir/merged.germline.indel --annotation-file $germline_dir/merged.germline.indel.transcript-annotation --output-file $germline_dir/merged.germline.indel.annotated";
		print SCRIPT "$cmd\n";
	}

	system("bsub -q long -R\"select[type==LINUX64 && mem>8000] rusage[mem=8000]\" -M 8000000 -oo $germline_dir/germline-indel.sh.log \"sh $germline_dir/germline-indel.sh\"");

}



################################################################################################
# Get SAMtools SNPs - get the SNPs from SAMtools filtered post-annotation file
#
################################################################################################

sub get_samtools_snps
{
	my ($variant_file, $output_file) = @_;	
	my $num_snps = 0;
	## Check for Tier 1 SNVs ##
	
	if(-e $variant_file)
	{
		open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";

		## Parse the Tier 1 SNVs file ##
	
		my $input = new FileHandle ($variant_file);
		my $lineCounter = 0;
	
		while (<$input>)
		{
			chomp;
			my $line = $_;
			$lineCounter++;

			my ($chrom, $chr_start, $chr_stop, $ref, $var, $var_type, $gene, $transcript, $organism, $source, $version, $strand, $status, $trv_type) = split(/\t/, $line);
			
			if($var_type eq "SNP" && is_tier1($trv_type))
			{
				$num_snps++;
				print OUTFILE "$line\n";
			}
		}
		
		close($input);
		
		close(OUTFILE);
	}
	
	return($num_snps);




}


###############################################################################################
# is_tier1 - return 1 if this trv_type is considered tier 1
#
################################################################################################

sub is_tier1
{
	my $trv_type = shift(@_);
	
	if($trv_type eq 'missense' || $trv_type eq 'nonsense' || $trv_type eq 'nonstop' || $trv_type eq 'splice_site' || $trv_type =~ 'frame')
	{
		return(1);
	}
	
	return(0);
}

################################################################################################
# Output Germline Files - output the files for loh analysis 
#
################################################################################################

sub output_loh_files
{
	my ($self, $build_dir, $loh_dir, $normal_bam, $tumor_bam) = @_;

	## Get VarScan Germline File ##
	
	my $germline_snv_file = `ls $build_dir/varScan.output.snp.formatted.Germline.hc`;
	chomp($germline_snv_file);
	my $germline_output_file = "$loh_dir/varScan.germline.snp";
	
	if($germline_snv_file && !(-e $germline_output_file))
	{	
		## Apply a preliminary filter ##	
		system("gmt varscan filter-variant-calls --variant $germline_snv_file --output $germline_output_file.unfiltered --min-normal-cov 8 --min-tumor-cov 8 --min-normal-freq 25 --max-normal-freq 75 --min-tumor-freq 0 --max-tumor-freq 100");
		
		## Apply the FP-filter ##
		my $cmd = "gmt somatic filter-false-positives --variant-file $germline_output_file.unfiltered --bam-file $tumor_bam --output-file $germline_output_file";
		$cmd .= " --reference " . $self->reference if($self->reference);
		system("bsub -q long -R\"select[type==LINUX64 && mem>8000 && tmp>2000] rusage[mem=8000]\" -M 8000000 $cmd");
	}
	

	## Get VarScan LOH File ##
	
	my $loh_snv_file = `ls $build_dir/varScan.output.snp.formatted.LOH.hc`;
	chomp($loh_snv_file);
	my $loh_output_file = "$loh_dir/varScan.loh.snp";
	
	if($loh_snv_file && !(-e $loh_output_file))
	{	
		## Apply a preliminary filter ##	
		system("gmt varscan filter-variant-calls --variant $loh_snv_file --output $loh_output_file.unfiltered --min-normal-cov 8 --min-tumor-cov 8 --min-normal-freq 25 --max-normal-freq 75 --min-tumor-freq 0 --max-tumor-freq 100");


		## Apply the FP-filter ##
		my $cmd = "gmt somatic filter-false-positives --variant-file $loh_output_file.unfiltered --bam-file $normal_bam --output-file $loh_output_file";
		$cmd .= " --reference " . $self->reference if($self->reference);
		system("bsub -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT} -R\"select[type==LINUX64 && mem>8000 && tmp>2000] rusage[mem=8000]\" -M 8000000 $cmd");
	}
	
	## If both files exist, process them ##
	my $merged_output_file = "$loh_dir/varScan.merged.snp";
	if(-e $germline_output_file && -e $loh_output_file)
	{
		print "Merging LOH and Germline calls...\n";
		system("cat $germline_output_file $loh_output_file >$merged_output_file");
		system("gmt capture sort-by-chr-pos --input $merged_output_file --output $merged_output_file");
		system("gmt varscan filter-variant-calls --variant $merged_output_file --output $merged_output_file.hc --min-normal-cov 8 --min-tumor-cov 8 --min-normal-freq 25 --max-normal-freq 75 --min-tumor-freq 0 --max-tumor-freq 100");
#		system("echo \"chrom\tpos\tnfreq\ttfreq\tstatus\" >$merged_output_file.hc.infile");
#		system("cut -f 1,2,8,12,14 $merged_output_file.hc | perl -pe 's/\%//g' >>$merged_output_file.hc.infile");
	}
}




################################################################################################
# Output Germline Files - output the files for loh analysis 
#
################################################################################################

sub process_loh
{
	my $loh_dir = shift(@_);

	my $germline_snp = "$loh_dir/varScan.germline.snp";
	my $loh_snp = "$loh_dir/varScan.loh.snp";
	
	if(-e $germline_snp && -e $loh_snp)
	{
		my $combined_snp = "$loh_dir/varScan.combined.snp";
		system("cat $germline_snp $loh_snp >$combined_snp");
		system("gmt capture sort-by-chr-pos --input $combined_snp --output $combined_snp");
		my $cmd = "gmt varscan loh-segments --variant-file $combined_snp --output-basename $combined_snp.loh";
		system("bsub -q short -R\"select[type==LINUX64 && mem>1000] rusage[mem=1000]\" \"$cmd\"");
	}
	else
	{
		warn "Warning: Germline/LOH SNP file(s) $germline_snp $loh_snp missing from $loh_dir\n";
	}
}


sub byChrPos
{
	my ($chr_a, $pos_a) = split(/\t/, $a);
	my ($chr_b, $pos_b) = split(/\t/, $b);
	
	$chr_a cmp $chr_b
	or
	$pos_a <=> $pos_b;
}


1;

