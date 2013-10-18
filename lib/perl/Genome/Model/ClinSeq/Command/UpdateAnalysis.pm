package Genome::Model::ClinSeq::Command::UpdateAnalysis;

#Written by Malachi Griffith

use strict;
use warnings;
use Genome;
use Genome::Model::ClinSeq;
use Data::Dumper;
use Term::ANSIColor qw(:constants);
use Time::Piece;

class Genome::Model::ClinSeq::Command::UpdateAnalysis {
    is => 'Command::V2',
    has_optional => [
        individual => { 
              is => 'Genome::Individual',
              is_many => 0,
              require_user_verify => 0,
              doc => 'Individual to query',
        },
        samples => {
              is => 'Genome::Sample',
              is_many => 1,
              require_user_verify => 0,
              doc => 'Sample(s) to target for clinseq analysis'
        },
        sample_type_filter => {
              is => 'Text',
              default => 'pcr product,pooled library,pooled dna,ipr product',
              doc => 'When displaying samples, filter out those with certain sample types. [comma separate list]',
        },
        _ref_align_pp_id => {
              is => 'Number',
              default => '2635769',
        },
        ref_align_pp => {
              is => 'Genome::ProcessingProfile',
              id_by => '_ref_align_pp_id',
              doc => 'Desired reference-alignment processing profile',
        },
        _wgs_somatic_variation_pp_id => {
              is => 'Number',
              default => '2762562',
        },
        wgs_somatic_variation_pp => {
              is => 'Genome::ProcessingProfile',
              id_by => '_wgs_somatic_variation_pp_id',
              doc => 'Desired WGS somatic-variation processing profile',
        },
        _exome_somatic_variation_pp_id => {
              is => 'Number',
              default => '2762563',
        },
        exome_somatic_variation_pp => {
              is => 'Genome::ProcessingProfile',
              id_by => '_exome_somatic_variation_pp_id',
              doc => 'Desired Exome somatic-variation processing profile',
        },
        _rnaseq_pp_id => {
              is => 'Number',
              default => '2762841',
        },
        rnaseq_pp => {
              is => 'Genome::ProcessingProfile',
              id_by => '_rnaseq_pp_id',
              doc => 'Desired rna-seq processing profile',
        },
        _differential_expression_pp => {
              is => 'Number',
              default => '2760181',
        },
        differential_expression_pp => {
              is => 'Genome::ProcessingProfile',
              id_by => '_differential_expression_pp',
              doc => 'Desired differential-expression processing profile',
        },
        _clinseq_pp_id => {
              is => 'Number',
              default => '2649924',
        },
        clinseq_pp => {
              is => 'Genome::ProcessingProfile',
              id_by => '_clinseq_pp_id',
              doc => 'Desired clin-seq processing profile',
        },
        _reference_sequence_build_id => {
              is => 'Number',
              default => '106942997',
        },
        reference_sequence_build => {
              is => 'Genome::Model::Build',
              id_by => '_reference_sequence_build_id',
              doc => 'Desired reference sequence build',
        },
        _annotation_build_id => {
              is => 'Number',
              default => '124434505',
        },
        annotation_build => {
              is => 'Genome::Model::Build',
              id_by => '_annotation_build_id',
              doc => 'Desired reference annotation build',
        },
        _dbsnp_build_id => {
              is => 'Number',
              default => '127786607',
        },
        dbsnp_build => {
              is => 'Genome::Model::Build',
              id_by => '_dbsnp_build_id',
              doc => 'Desired dbSNP build',
        },
        _previously_discovered_variations_id => {
              is => 'Number',
              default => '127786607',
        },
        previously_discovered_variations => {
              is => 'Genome::Model::Build',
              id_by => '_previously_discovered_variations_id',
              doc => 'Desired previously discovered variants build',
        },
        display_defaults => {
              is => 'Boolean',
              doc => 'Display current default processing profiles and annotation/reference genome inputs',
        },
        display_samples => {
              is => 'Boolean',
              doc => 'Display sample info and abort',
        },
        normal_sample_common_names => {
              #TODO: Is there a better way to determine which samples are 'normal'?
              is => 'Text',
              default => 'normal',
              doc => 'The possible sample common names used in the database to specify a Normal sample',
        },
        tumor_sample_common_names => {
              #TODO: Is there a better way to determine which samples are 'tumor'?
              is => 'Text',
              default => 'tumor|met|post treatment|recurrence met|pre-treatment met|pin lesion|relapse|xenograft',
              doc => 'The possible sample common names used in the database to specify a Tumor sample',
        },
        instrument_data_to_exclude => {
              is => 'Text',
              doc => 'Instrument data to exclude from all consideration. Supply as a comma separated list of instrument data IDs. Supply a mix of dna and rna data if required.',
        },
        skip_check_archived => {
              is => 'Boolean',
              doc => 'Check if builds are currently archived',
        },
        check_archivable_status => {
              is => 'Boolean',
              doc => 'Warn if builds are currently set to be *archivable*.  Even if false you will still be warned if something is already archived',
        },
        force => {
              is => 'Boolean',
              doc => 'Allow certain warnings/errors to be by-passed',
        },
        ignore_models_matching => {
              is => 'Text',
              is_optional => 1,
              doc => "expression matching models which should be ignored",
        },
   ],
    doc => 'evaluate models/builds for an individual and help create/update a clinseq model that meets requested criteria',
};

sub help_synopsis {
    return <<EOS

Creating or updating a clin-seq analysis involves the following three basic steps:

Step 1. Summarize the current default processing profiles and inputs:
genome model clin-seq update-analysis  --display-defaults

Step 2. Examine available samples for your individual:
genome model clin-seq update-analysis  --individual='2878747495'
genome model clin-seq update-analysis  --individual='H_KA-306905'
genome model clin-seq update-analysis  --individual='common_name=AML103'

Step 3: Run the update-analysis component to determine what models need to be created:
genome model clin-seq update-analysis  --individual='H_KA-306905' --samples='id in [2878747496,2878747497,2879495575]'
genome model clin-seq update-analysis  --individual='H_KA-306905' --samples='name in ["H_KA-306905-1121472","H_KA-306905-1121474","H_KA-306905-S.4294"]'

All processing-profile and build input parameters can be specified by ID, name, etc. and should resolve

EOS
}

sub help_detail {
    return <<EOS
For a given individual, find the available samples and display to the user to select those desired for clinseq analysis

A clin-seq model can consist of various combinations of reference alignment, wgs somatic, exome somatic, and rna-seq models, 

Only up to two DNA and two RNA samples may be specified. You can run with only DNA, only RNA or both.

Once samples are set by the user, search for reference alignment, somatic variation, and rna-seq models for these samples

Determine whether these models meet certain criteria with respect to inputs and processing profile selections

EOS
}


sub execute {
  my $self = shift;

  #If the user selected the --display-defaults option, simply print out a summmary and exit
  if ($self->display_defaults){
    $self->status_message("\n\nSummarizing default processing profiles and inputs");
    $self->display_processing_profiles;
    $self->display_inputs;
    $self->status_message("\nExiting...  --display-defaults mode is used for summarizing purposes only\n\n");
    return 1;
  }

  #Make sure an individual is defined
  unless ($self->individual){
    $self->error_message("Missing required parameter: --individual.");
    exit 1;
  }

  my $individual = $self->individual;

  #Display basic info about the individual
  $self->display_individual;

  #Summarize the desired processing profiles (either supplied by user, or otherwise defaults): reference-alignment, somatic-variation (WGS & Exome), rna-seq, clin-seq
  $self->display_processing_profiles;

  #Summarize the desired input parameters (either supplied by user, or otherwise defaults): reference sequence build, reference annotation build, etc.
  $self->display_inputs;

  #Get samples
  my @samples;
  if ($self->samples){
    #If the user supplied a list of target sample IDs or names, retrieve the sample objects
    @samples = $self->get_samples;
  }else{
    #If the user did not display the samples to use for this individual, display those available and exit.
    $self->get_samples;
    $self->status_message("\nUser must select which samples with --samples in order to proceed\n\n");
    return 1;
  }

  #If the user selected the --display-samples option, simply print out a summmary and exit
  if ($self->display_samples){
    $self->status_message("\nExiting...  --display-samples mode is used for summarizing purposes only\n\n");
    return 1;
  }

  #Get the subset of samples that are of the type DNA or RNA
  $self->status_message("\nGET SAMPLES BY TYPE");
  my @dna_samples = $self->dna_samples('-samples'=>\@samples);
  my @rna_samples = $self->rna_samples('-samples'=>\@samples);

  #Once the user has supplied the desired samples, proceed with additional examination of instrument data, models, etc.
  #In the following checks, a 'suitable' model uses the desired processing profile and desired input parameters
  #Each method will return a status indicating whether it is safe to proceed to the next stage
  #It is safe to proceed if:
  #  - 1. the neccessary data does not even exist
  #  - 2. the model already exists and meets all the desired criteria
  #It is not safe to proceed if:
  #  - 1. a new model needs to be created (in this case, do this first)
  #  - 2. a model exists that meets the criteria but it has not succeeded yet
  $self->status_message("\nEXAMINE MODELS");

  #TODO: Are there suitable genotype microarray models in existence (If not, create)?  If so, what is the status?
  #TODO: Refer to this subroutine for guidance: get_genotype_microarray_model_id
  #TODO: Rework this logic to be more similar to the following steps.  Check for microarray models, create if they don't exist, return all matching models if they do
  #TODO: Then feed these as inputs into the next step that checks ref-align models.  If there are multiple matching models get the 'best' one.
  #TODO: If builds are underway, don't proceed to the next step.  If there is no microarray data proceed anyway.
  #TODO: See docs for details: genome model define genotype-microarray --help


  #Are there suitable WGS reference alignment models in existence (If not, create)?  If so, what is the status?
  #- Are there tumor/normal DNA samples?
  #- Is there tumor/normal WGS instrument data? If so, is it all included in the existing model?
  my @normal_wgs_ref_align_models;
  my @tumor_wgs_ref_align_models;
  if (scalar(@dna_samples) >= 1){
    $self->status_message("\nWGS REFERENCE-ALIGNMENT MODELS");
    @normal_wgs_ref_align_models = $self->check_ref_align_models('-data_type'=>'wgs', '-tissue_type'=>$self->normal_sample_common_names, '-dna_samples'=>\@dna_samples);
    @tumor_wgs_ref_align_models = $self->check_ref_align_models('-data_type'=>'wgs', '-tissue_type'=>$self->tumor_sample_common_names, '-dna_samples'=>\@dna_samples);
  }

  my @wgs_somatic_variation_models;
  if (scalar(@dna_samples) == 2){
    #Is there a suitable WGS somatic variation model in existence (If not, create)?  If so, what is the status?
    #- Only proceed with this if the prerequisite WGS tumor/normal ref-align models exist
    if (scalar(@normal_wgs_ref_align_models) && scalar(@tumor_wgs_ref_align_models)){
      $self->status_message("\nWGS SOMATIC-VARIATION MODELS");
      @wgs_somatic_variation_models = $self->check_somatic_variation_models('-data_type'=>'wgs', '-dna_samples'=>\@dna_samples, '-normal_models'=>\@normal_wgs_ref_align_models, '-tumor_models'=>\@tumor_wgs_ref_align_models);
    }
  }

  #Are there suitable Exome reference alignment models in existence (If not, create)? If so, what is the status?
  #- Are there tumor/normal DNA samples?
  #- Is there tumor/normal Exome instrument data?
  #- What is the target region set name (TRSN) and region of interest (ROI)? If there is data, is it all included in the existing model?
  $self->status_message("\nEXOME REFERENCE-ALIGNMENT MODELS");
  my @normal_exome_ref_align_models;
  my @tumor_exome_ref_align_models;
  if (scalar(@dna_samples) >= 1){
    @normal_exome_ref_align_models = $self->check_ref_align_models('-data_type'=>'exome', '-tissue_type'=>$self->normal_sample_common_names, '-dna_samples'=>\@dna_samples);
    @tumor_exome_ref_align_models = $self->check_ref_align_models('-data_type'=>'exome', '-tissue_type'=>$self->tumor_sample_common_names, '-dna_samples'=>\@dna_samples);
  }

  my @exome_somatic_variation_models;
  if (scalar(@dna_samples) == 2){
    #Is there a suitable Exome somatic variation model in existence
    #- Only proceed with this if the prerequisite Exome tumor/normal ref-align models exist
    if (scalar(@normal_exome_ref_align_models) && scalar(@tumor_exome_ref_align_models)){
      $self->status_message("\nEXOME SOMATIC-VARIATION MODELS");
      @exome_somatic_variation_models = $self->check_somatic_variation_models('-data_type'=>'exome', '-dna_samples'=>\@dna_samples, '-normal_models'=>\@normal_exome_ref_align_models, '-tumor_models'=>\@tumor_exome_ref_align_models);
    }
  }
  
  #Are there suitable rna-seq models in existence (If not, create)?  If so, what is the status?  
  #- Are there tumor/normal RNA samples?
  #- Is there tumor/normal RNA-seq instrument data? If so, is it all included in the existing model?
  my @normal_rnaseq_models;
  my @tumor_rnaseq_models;
  if (scalar(@rna_samples)){
    $self->status_message("\nRNA-SEQ MODELS");
    @normal_rnaseq_models = $self->check_rnaseq_models('-tissue_type'=>$self->normal_sample_common_names, '-rna_samples'=>\@rna_samples);
    @tumor_rnaseq_models = $self->check_rnaseq_models('-tissue_type'=>$self->tumor_sample_common_names, '-rna_samples'=>\@rna_samples);
  }

  #Gather the 'best' suitable WGS somatic-variation, Exome somatic-variation, tumor RNA-seq, and normal RNA-seq models
  my $best_wgs_somatic_variation_model = $self->select_best_model('-models'=>\@wgs_somatic_variation_models);
  my $best_exome_somatic_variation_model = $self->select_best_model('-models'=>\@exome_somatic_variation_models);
  my $best_normal_rnaseq_model = $self->select_best_model('-models'=>\@normal_rnaseq_models);
  my $best_tumor_rnaseq_model = $self->select_best_model('-models'=>\@tumor_rnaseq_models);

  #If there is a 'best' tumor AND normal rna-seq model, check if there is a differential-expression model with these two
  my @de_models;
  if ($best_normal_rnaseq_model && $best_tumor_rnaseq_model){
    $self->status_message("\nDIFFERENTIAL-EXPRESSION MODELS");
    @de_models = $self->check_de_models('-normal_rnaseq_model'=>$best_normal_rnaseq_model, '-tumor_rnaseq_model'=>$best_tumor_rnaseq_model);
  }
  my $best_de_model = $self->select_best_model('-models'=>\@de_models);

  #Is there a suitable clin-seq model in existence that uses these models (If not, create it)?  
  #If there is a matching clin-seq model, what is the status?  Remind the user whether this clinseq model is currently set as do-not-archive...
  $self->status_message("\nCLIN-SEQ MODELS");
  my $clinseq_model = $self->check_clinseq_models('-wgs_somatic_variation_model'=>$best_wgs_somatic_variation_model,
                                                  '-exome_somatic_variation_model'=>$best_exome_somatic_variation_model,
                                                  '-normal_rnaseq_model'=>$best_normal_rnaseq_model,
                                                  '-tumor_rnaseq_model'=>$best_tumor_rnaseq_model,
                                                  '-de_model'=>$best_de_model,
                                                  '-samples'=>\@samples);


  if ($clinseq_model){
    my $clinseq_inputs_ok = $self->check_clinseq_inputs('-model'=>$clinseq_model);
    if ($clinseq_inputs_ok){
      my $clinseq_model_id = $clinseq_model->id;
      my $clinseq_model_name = $clinseq_model->name;
      my @inputs;
      push (@inputs, 'wgs') if $clinseq_model->wgs_model;
      push (@inputs, 'exome') if $clinseq_model->exome_model;
      push (@inputs, 'tumor_rnaseq') if $clinseq_model->tumor_rnaseq_model;
      push (@inputs, 'normal_rnaseq') if $clinseq_model->normal_rnaseq_model;
      push (@inputs, 'differential_expression') if $clinseq_model->de_model;
      my $input_string = join (", ", @inputs);
      $self->status_message("\nRequested clin-seq model exists and is complete:\n\tWith inputs: $input_string\n\t$clinseq_model_name <$clinseq_model_id>");
      my $date_completed = $clinseq_model->last_complete_build->date_completed;
      $date_completed =~ s/\.\d+$//;
      $self->status_message("\nThe last completed build finished running at " . $date_completed . ". If this was a while ago you might want to rebuild to take advantage of analysis improvements...\n\n");
    }
  }else{
    $self->status_message("\nRequested clin-seq model is NOT ready - see above for instructions\n\n");
  }

  return 1;
}


#Resolve a human readable individual name from an individual object
sub get_final_individual_name{
  my $self = shift;
  my $individual = $self->individual;
  my $final_name;
  my $upn = $individual->upn;
  if ($individual->common_name){
    $final_name = $individual->common_name;
  }elsif($individual->name){
    $final_name = $individual->name;
  }else{
    $final_name = $individual->id;
  }
  if ($upn){
    $final_name .= " ($upn)";
  }
  return $final_name;
}


#Print a brief display of details pertaining to the individual
sub display_individual{
  my $self = shift;
  my $individual = $self->individual;
  my $id = $individual->id;
  my $species_name =  $individual->species_name;
  my $name = $individual->name || "NULL";
  my $upn = $individual->upn || "NULL";
  my $common_name = $individual->common_name || "NULL";
  my $gender = $individual->gender || "NULL";
  my $ethnicity = $individual->ethnicity || "NULL";
  my $race = $individual->race || "NULL";
  $self->status_message("\nEXAMINE INDIVIDUAL");
  $self->status_message("id\tname\tcommon_name\tspecies_name\tupn\tgender\tethnicity\trace");
  $self->status_message("$id\t$name\t$common_name\t$species_name\t$upn\t$gender\t$ethnicity\t$race");
  return;
}


#Get samples for the individual and filter out those with the specified sample type
sub get_samples{
  my $self = shift;
  my $individual = $self->individual;

  my @samples;
  my @final_samples;
  if ($self->samples){
    @samples = $self->samples;
  }else{
    @samples = $individual->samples;
  }
  @samples = sort { $a->name cmp $b->name } @samples;
  
  my $sample_count = scalar(@samples);
  $self->status_message("\nEXAMINE SAMPLES");
  $self->status_message("Found $sample_count samples:");
  $self->status_message("id\tname\tsample_common_name\tsample_type\tcell_type\ttissue_desc\tdefault_genotype_data_id\tmodel_count\tlibrary_count\tid_count");
  my $skip_count = 0;
  my $sample_mismatch = 0;
  foreach my $sample (@samples){
    my $id = $sample->id;
    my $name = $sample->name;
    my $sample_common_name = $sample->common_name || "NULL";
    my $extraction_label = $sample->extraction_label || "NULL";
    my $extraction_type = $sample->extraction_type || "NULL";
    my $sample_type = $sample->sample_type || "NULL";
    my $extraction_desc = $sample->extraction_desc || "NULL";
    my $cell_type = $sample->cell_type || "NULL";
    my $tissue_label = $sample->tissue_label || "NULL";
    my $tissue_desc = $sample->tissue_desc || "NULL";
    my $organ_name = $sample->organ_name || "NULL";
    my $disease = $sample->disease || "NULL";
    my $default_genotype_data_id = $sample->default_genotype_data_id || "NULL";
    my @models = $sample->models;
    my $model_count = scalar(@models);
    my @libraries = $sample->libraries;
    my $library_count = scalar(@libraries);
    my @instrument_data = $sample->instrument_data;
    @instrument_data = @{$self->exclude_instrument_data('-instrument_data'=>\@instrument_data)};    
    my $id_count = scalar(@instrument_data);
    my $skip = 0;
    if ($self->sample_type_filter){
      my @stf = split(",", $self->sample_type_filter);
      foreach my $stf (@stf){
        if ($stf =~ /$sample_type/i){
          $skip_count++;
          $skip = 1;
        }
      }
    }
    next if $skip;

    #Make sure the individual specified matches that of every sample
    my $patient_id = $sample->patient->id;
    my $individual_id = $individual->id;
    unless ($patient_id eq $individual_id){
      $self->error_message("ID of individual supplied by user ($individual_id) does not match that associated with a sample ($patient_id)");
      $sample_mismatch++;
    }
    #$self->status_message("id\tname\tsample_common_name\tsample_type\tcell_type\ttissue_desc\tdefault_genotype_data_id\tmodel_count\tlibrary_count\tid_count");
    $self->status_message("$id\t$name\t$sample_common_name\t$sample_type\t$cell_type\t$tissue_desc\t$default_genotype_data_id\t$model_count\t$library_count\t$id_count");
    push(@final_samples, $sample);
  }

  if ($sample_mismatch){
    if ($self->force){
      $self->warning_message("Found $sample_mismatch samples provided by the user that do not match the specified patient.  Allowing since --force was used\n");
    }else{
      $self->warning_message("Found $sample_mismatch samples provided by the user that do not match the specified patient.  Aborting ...\n");
      exit 1;
    }
  }

  if ($self->sample_type_filter && $skip_count){
    $self->status_message("\nSkipped $skip_count samples due to matches to " . $self->sample_type_filter);
  }
  return @final_samples;
}


#Display the processing profiles to be targeted
sub display_processing_profiles{
  my $self = shift;
  $self->status_message("\nEXAMINE DESIRED PROCESSING PROFILES");
  $self->status_message("reference-alignment: " . $self->ref_align_pp->name . " <" . $self->ref_align_pp->id . ">");
  $self->status_message("wgs somatic-variation: " . $self->wgs_somatic_variation_pp->name . " <" . $self->wgs_somatic_variation_pp->id . ">");
  $self->status_message("exome somatic-variation: " . $self->exome_somatic_variation_pp->name . " <" . $self->exome_somatic_variation_pp->id . ">");
  $self->status_message("rna-seq: " . $self->rnaseq_pp->name . " <" . $self->rnaseq_pp->id . ">");
  $self->status_message("differential-expression: " . $self->differential_expression_pp->name . " <" . $self->differential_expression_pp->id . ">");
  $self->status_message("clin-seq: " . $self->clinseq_pp->name . " <" . $self->clinseq_pp->id . ">");
  return;
}


#Display the model inputs to be targeted
sub display_inputs{
  my $self = shift;
  $self->status_message("\nEXAMINE DESIRED MODEL INPUTS");
  $self->status_message("reference_sequence_build: " . $self->reference_sequence_build->__display_name__);
  $self->status_message("annotation_build: " . $self->annotation_build->name . " (" . $self->annotation_build->id . ")");
  $self->status_message("dbsnp_build: " . $self->dbsnp_build->__display_name__ . " (version " . $self->dbsnp_build->version . ")");
  $self->status_message("previously_discovered_variations: " . $self->previously_discovered_variations->__display_name__ . " (version " . $self->previously_discovered_variations->version . ")");

  my $cancer_annotation_db_id = Genome::Model::ClinSeq->__meta__->property("cancer_annotation_db")->default_value;
  my $cancer_annotation_db = Genome::Db::Tgi::CancerAnnotation->get($cancer_annotation_db_id);
  $self->status_message("cancer_annotation_db: " . $cancer_annotation_db_id . " (" . $cancer_annotation_db->data_directory . ")");

  my $misc_annotation_db_id = Genome::Model::ClinSeq->__meta__->property("misc_annotation_db")->default_value;
  my $misc_annotation_db = Genome::Db::Tgi::MiscAnnotation->get($misc_annotation_db_id);
  $self->status_message("misc_annotation_db: " . $misc_annotation_db_id . " (" . $misc_annotation_db->data_directory . ")");

  my $cosmic_annotation_db_id = Genome::Model::ClinSeq->__meta__->property("cosmic_annotation_db")->default_value;
  my $cosmic_annotation_db = Genome::Db::Cosmic->get($cosmic_annotation_db_id);
  $self->status_message("cosmic_annotation_db: " . $cosmic_annotation_db_id . " (" . $cosmic_annotation_db->data_directory . ")");

  #Make sure none of the basic input models/builds have been archived before proceeding...
  unless ($self->skip_check_archived){
    if ($self->reference_sequence_build->is_archived){
      $self->error_message("Reference sequencing build " . $self->reference_sequence_build->__display_name__ . " has been archived!");
      exit 1;
    }
    if ($self->annotation_build->is_archived){
      $self->error_message("Annotation build " . $self->annotation_build->name . " has been archived!");
      exit 1;
    }
    if ($self->dbsnp_build->is_archived){
      $self->error_message("dbSNP build " . $self->dbsnp_build->__display_name__ . " has been archived!");
      exit 1;
    }
  }

  return;
}


#Test for existence of DNA samples and return array of sample objects if found
sub dna_samples{
  my $self = shift;
  my %args = @_;
  my @samples = @{$args{'-samples'}};
    
  my $normal_sample_common_names = $self->normal_sample_common_names;
  my $tumor_sample_common_names = $self->tumor_sample_common_names;
  my @dna_samples;
  foreach my $sample (@samples){
    if ($sample->sample_type =~ /dna/i){
      push(@dna_samples, $sample);
    }
  }
  my $dna_sample_count = scalar(@dna_samples);
  $self->status_message("Found " . $dna_sample_count . " DNA samples");

  #Clin-Seq does not support more than 2 DNA samples and expects that one will be 'tumor' and one will be 'normal'
  if ($dna_sample_count > 2){
    $self->error_message("ClinSeq does not support more than 2 DNA samples ... please supply up to 1 tumor and 1 normal");
    exit 1;
  }
  my @normal_samples;
  my @tumor_samples;

  foreach my $s (@dna_samples){
    my $current_scn = $s->common_name || "NULL";
    push (@normal_samples, $s) if ($current_scn =~ /$normal_sample_common_names/i);
    push (@tumor_samples, $s) if ($current_scn =~ /$tumor_sample_common_names/i);
  }
  if (scalar(@normal_samples) > 1){
    $self->error_message("More than one normal DNA sample was specified for this individual - check samples or normal/tumor definitions");
  }
  if (scalar(@tumor_samples) > 1){
    $self->error_message("More than one tumor DNA sample was specified for this individual - check samples or normal/tumor definitions");
  }

  return (@dna_samples);
}


#Check for existence of RNA samples and return array of sample objects if found
sub rna_samples{
  my $self = shift;
  my %args = @_;
  my @samples = @{$args{'-samples'}};

  my $normal_sample_common_names = $self->normal_sample_common_names;
  my $tumor_sample_common_names = $self->tumor_sample_common_names;
  my @rna_samples;
  foreach my $sample (@samples){
    if ($sample->sample_type =~ /rna/i){
      push(@rna_samples, $sample);
    }elsif($sample->is_rna){
      push(@rna_samples, $sample);
    }
  }
  my $rna_sample_count = scalar(@rna_samples);
  $self->status_message("Found " . $rna_sample_count . " RNA samples");
  
  #Clin-Seq does not support more than 2 RNA samples and expects that one will be 'tumor' and one will be 'normal'
  if ($rna_sample_count > 2){
    $self->error_message("ClinSeq does not support more than 2 RNA samples ... please supply up to 1 tumor and 1 normal");
    exit 1;
  }
  my @normal_samples;
  my @tumor_samples;

  foreach my $s (@rna_samples){
    my $current_scn = $s->common_name || "NULL";
    push (@normal_samples, $s) if ($current_scn =~ /$normal_sample_common_names/i);
    push (@tumor_samples, $s) if ($current_scn =~ /$tumor_sample_common_names/i);
  }
  if (scalar(@normal_samples) > 1){
    $self->error_message("More than one normal RNA sample was specified for this individual - check samples or normal/tumor definitions");
  }
  if (scalar(@tumor_samples) > 1){
    $self->error_message("More than one tumor RNA sample was specified for this individual - check samples or normal/tumor definitions");
  }

  return (@rna_samples);
}


#Get all instrument data for the sample of the specified type (wgs or exome)
sub get_dna_instrument_data{
  my $self = shift;
  my %args = @_;
  my $sample = $args{'-sample'};
  my $data_type = $args{'-data_type'};

  my @sample_instrument_data = $sample->instrument_data;
  @sample_instrument_data = @{$self->exclude_instrument_data('-instrument_data'=>\@sample_instrument_data)};
  my $instrument_data_count = scalar(@sample_instrument_data);
  
  my @exome;
  my @wgs;
  my @unknown;
  my @other;
  my %trsns;

  foreach my $instrument_data (@sample_instrument_data){
    my $trsn = $instrument_data->target_region_set_name;
    if ($trsn){
      my $fl = Genome::FeatureList->get(name => $trsn);
      if (not $fl or not $fl->content_type) {
        push @unknown, $instrument_data;
      }elsif ($fl->content_type eq 'exome') {
        push @exome, $instrument_data;
        $trsns{$trsn}=1;
      }else {
        push @other, $instrument_data;
      }
    }else{
      push @wgs, $instrument_data;
    }
  }

  if ($data_type eq 'wgs'){
    $self->status_message("\tCould not find any wgs data") unless (scalar(@wgs));
    return @wgs;
  }elsif($data_type eq 'exome'){
    $self->status_message("\tCould not find any exome data") unless (scalar(@exome));
    return @exome;
  }elsif($data_type eq 'validation'){
    $self->status_message("\tCould not find any exome data") unless (scalar(@exome));
    return @exome;
  }else{
    $self->error_message("Data type not understood in get_dna_instrument_data");
  }
}


#Determine exome vs WGS of a model.  WGS models do not have exome data, Exome models have at least one lane of exome data
    
sub determine_model_data_type{
  my $self = shift;
  my %args = @_;
  my $model = $args{'-model'};
  my $model_data_type = "unknown";
  my @model_instrument_data = $model->instrument_data;
  my $instrument_data_count = scalar(@model_instrument_data);

  my $wgs_lane_count = 0;
  my $exome_lane_count = 0;
  my $other_lane_count = 0;
  my @unknown;
  my @other;
  foreach my $instrument_data (@model_instrument_data){
    my $trsn = $instrument_data->target_region_set_name;
    if ($trsn){
      my $fl = Genome::FeatureList->get(name => $trsn);
      if (not $fl or not $fl->content_type) {
        push @unknown, $instrument_data;
      }elsif ($fl->content_type eq 'exome') {
        $exome_lane_count++;
      }else {
        push @other, $instrument_data;
        $other_lane_count++;
      }
    }else{
      $wgs_lane_count++;
    }
  }

  if ($instrument_data_count == $wgs_lane_count){
    $model_data_type = "wgs";
  }elsif($exome_lane_count >= 1){
    $model_data_type = "exome";
  } 

  return $model_data_type;
}


#Determine the target region set name of the instrument data associated with the model and make sure the same values were specified in creating the model for both TRSN and ROI
sub check_model_trsn_and_roi{
  my $self = shift;
  my %args = @_;
  my $model = $args{'-model'};
  my $data_type = $args{'-data_type'};

  #If this not exome data, automatically pass the model
  return 1 unless ($data_type eq 'exome');

  my %trsns;
  my @model_instrument_data = $model->instrument_data;
  my $model_trsn = $model->target_region_set_name;
  my $model_roi = $model->region_of_interest_set_name;
  my $model_id = $model->id; 

  my $trsn_ref;
  foreach my $instrument_data (@model_instrument_data){
    my $trsn = $instrument_data->target_region_set_name;
    if ($trsn){
      $trsns{$trsn}=1;
      $trsn_ref = $trsn;
    }
  }
  
  #Watch out for cases where multiple TRSNs have been combined...
  my $trsn_count = keys %trsns;
  if ($trsn_count >= 2){
    $self->warning_message("Intrument data from more than one target region set are being combined... (model id = $model_id)");
  }elsif($trsn_count == 0){
    $self->error_message("There is no instrument data with a target region set name!  How is this an exome data set?");
    exit 1;
  }

  #The target region set name of the model should match that of the data.
  my $trsn_match = 0;
  if ($trsn_ref eq $model_trsn){
    $trsn_match = 1;
  }

  #The ROI will normally match but may not if the capture reagent was designed on a previous version of the human genome...
  #Not sure how to deal with this in a good way that is automatic...
  #In the simplest case the desired region_of_interest_name is the same as the desired target_region_set_name
  #In some specific case, certain differences are desired depending on the reference genome version being used for example
  #Hard code these exceptions here
  my $roi_ref = $self->get_roi_name('-target_region_set_name'=>$trsn_ref);

  my $roi_match = 0;
  if ($roi_ref eq $model_roi){
    $roi_match = 1;
  }

  if ($trsn_match && $roi_match){
    return 1;
  }else{
    return undef;
  }
}


#For an array of instrument data objects, determine the target_region_set_name and perform basic checks
sub get_trsn{
  my $self = shift;
  my %args = @_;
  my @instrument_data = @{$args{'-instrument_data'}};

  my %trsns;
  my $trsn_ref;

  foreach my $instrument_data (@instrument_data){
    my $trsn = $instrument_data->target_region_set_name;
    if ($trsn){
      $trsns{$trsn}=1;
      $trsn_ref = $trsn;
    }
  }
  
  #Watch out for cases where multiple TRSNs have been combined...
  my $trsn_count = keys %trsns;
  if ($trsn_count >= 2){
    $self->warning_message("Intrument data from more than one target region set are being combined...");
  }elsif($trsn_count == 0){
    $self->error_message("There is no instrument data with a target region set name!  How is this an exome data set?");
    exit 1;
  }
  return $trsn_ref; 
}


#Determine the desired region of interest name based on target region set name and reference alignment version
sub get_roi_name{
  my $self = shift;
  my %args = @_;
  my $trsn = $args{'-target_region_set_name'};

  my $roi_name = $trsn;
  if ($self->reference_sequence_build->name =~ /GRCh37/i && $trsn eq 'hg18 nimblegen exome version 2'){    
    $roi_name = "hg19 nimblegen exome version 2";
  }elsif ($self->reference_sequence_build->name =~ /GRCh37/i && $trsn eq 'SeqCap EZ Human Exome v2.0'){
     $roi_name = "hg19 nimblegen exome version 2";
  }
  return $roi_name;
}


#Compare instrument data on a model to that available for the sample (exome or wgs) and warn if data is missing
sub check_for_missing_and_excluded_data{
  my $self = shift;
  my %args = @_;
  my @models = @{$args{'-models'}};
  my @sample_instrument_data = @{$args{'-sample_instrument_data'}};

  #Before doing anything, exclude models that have any instrument data that is supposed to be excluded
  if ($self->instrument_data_to_exclude){
    my @models_without_excluded_data;
    
    #Check format of list of data to exclude
    my @exclude_list = split(",", $self->instrument_data_to_exclude);
    unless (scalar(@exclude_list)){
      $self->error_message("Could not obtain instrument data IDs from list supplied by --instrument_data_to_exclude");
      exit 1;
    }

    foreach my $model (@models){
      my $found_excluded_data = 0;
      my @model_instrument_data = $model->instrument_data;

      foreach my $model_instrument_data (@model_instrument_data){
        my $mid = $model_instrument_data->id;
        foreach my $eid (@exclude_list){
          if ($eid eq $mid){
            $found_excluded_data++;
          }
        }
      }
      unless ($found_excluded_data){
        push(@models_without_excluded_data, $model);
      }
    }
    @models = @models_without_excluded_data;
  }

  my @final_models1;
  my $complete_model = 0; #Is there at least one model that is not missing any data?
  my %problem_models;
  foreach my $model (@models){
    my @model_instrument_data = $model->instrument_data;
    my @missing_model_data;

    foreach my $sample_instrument_data (@sample_instrument_data){
      my $sid = $sample_instrument_data->id;
      my $match = 0;
      foreach my $model_instrument_data (@model_instrument_data){
        my $mid = $model_instrument_data->id;
        $match = 1 if ($mid eq $sid);
      }
      push(@missing_model_data, $sid) unless $match;
    }
    if (scalar(@missing_model_data)){
      my $model_id = $model->id;
      my $id_string = join(",", @missing_model_data);
      $problem_models{$model_id}{id_string} = $id_string;
      $problem_models{$model_id}{name} = $model->name;
    }else{
      push(@final_models1, $model);
    }
  }

  #If no suitable models were found - report on those missing instrument data
  unless (scalar(@final_models1)){
    foreach my $model_id (sort keys %problem_models){
      my $id_string = $problem_models{$model_id}{id_string};
      my $model_name = $problem_models{$model_id}{name};
      $self->status_message("\tName: $model_name ($model_id)");  
      $self->status_message("\tWARNING -> Model: $model_id appears to be missing the following instrument data: $id_string");
      $self->status_message("\tYou should consider performing the following update before proceeding:");
      $self->status_message("\tgenome model instrument-data assign --instrument-data='$id_string'  --model=$model_id\n\t\t\tgenome model build start $model_id");
    }
  }

  #Does the lastest build actually use all the data?
  my @final_models2;
  my $complete_build = 0; #Is there at least one model where the last build is not missing any data?
  my %problem_builds;
  foreach my $model (@final_models1){
    my $last_build = $model->latest_build;
    if ($last_build){
      my @build_instrument_data = $last_build->instrument_data;
      my @missing_build_data;
      foreach my $sample_instrument_data (@sample_instrument_data){
        my $sid = $sample_instrument_data->id;
        my $match = 0;
        foreach my $build_instrument_data (@build_instrument_data){
          my $mid = $build_instrument_data->id;
          $match = 1 if ($mid eq $sid);
        }
        push(@missing_build_data, $sid) unless $match;
      }
      if (scalar(@missing_build_data)){
        my $model_id = $model->id;
        my $id_string = join(",", @missing_build_data);
        $problem_builds{$model_id}{id_string} = $id_string;
        $problem_builds{$model_id}{name} = $model->name;
      }else{
        push(@final_models2, $model);
      }
    }
  }

  #Warn about builds missing intrument data
  unless (scalar(@final_models2)){
    foreach my $model_id (sort keys %problem_builds){
      my $id_string = $problem_builds{$model_id}{id_string};
      $self->status_message("\tWARNING -> Last build of Model: $model_id appears to be missing the following instrument data: $id_string");
      $self->status_message("\tIt is on the model so you may need to start a new build before proceeding:");
      $self->status_message("\tgenome model build start $model_id");
    }
  }

  #Warn about builds that have all the data but are stalled on a non-succeeded status
  foreach my $model (@final_models2){
    my $model_id = $model->id;
    my $latest_build = $model->latest_build;
    unless ($latest_build->status eq "Succeeded"){
      $self->status_message("\tWARNING -> Last build of Model: $model_id has a status of: " . $latest_build->status);
    }
  }

  return (\@final_models2);
}


#Obtain a genotype microarray model object (if it exists) for a particular sample
sub get_genotype_microarray_model_id{
  my $self = shift;
  my %args = @_;
  my $sample = $args{'-sample'};

  my $genotype_microarray_model_id = 0;

  #Get the imported data that should have been used for the genotype microarray model
  my $default_genotype_data = $sample->default_genotype_data;

  #If there is no default genotype data defined, return 0
  return $genotype_microarray_model_id unless $default_genotype_data;

  #If there is genotype microarray data, look for GenotypeMicroarray models
  my @models = $sample->models;
  my @final_models;
  my @skipped;
  foreach my $model (@models){
    next unless ($model->class eq "Genome::Model::GenotypeMicroarray");

    unless ($model->reference_sequence_build) {
        push @skipped, "No reference_sequence_build on microarray model " . $model->__display_name__;
        next;
    }
    my $microarray_model_id = $model->id;

    #Make sure the genotype microarray model is on the specified version of the reference genome
    next unless ($model->reference_sequence_build->id eq $self->reference_sequence_build->id);

    #TODO: Make sure the genotype microarray model is using the specified version of dbSNP?  Skip this test for now...
    #if ($model->can("dbsnp_build")){
    #  next unless ($model->dbsnp_build->id eq $self->dbsnp_build->id);
    #}else{
    #  next;
    #}
    push (@final_models, $model);
  }

  # only omit warnings if zero things were found
  if (@final_models == 0) {
    if (@skipped) {
      for my $msg (@skipped) {
        $self->warning_message($msg)
      }
    }
  }

  if (scalar(@final_models)){
    #If there are multiple suitable genotype microarray models, select the best one
    my $best_model = $self->select_best_model('-models'=> \@final_models);
    $genotype_microarray_model_id = $best_model->id;
    return $genotype_microarray_model_id;
  }else{
    return 0;
  }
}


#Gather reference alignment models for a single subject that meet all the specified criteria
sub check_ref_align_models{
  my $self = shift;
  my %args = @_;
  my $data_type = $args{'-data_type'};
  my $tissue_type = $args{'-tissue_type'};
  my @dna_samples = @{$args{'-dna_samples'}};

  my @final_models;

  #Make sure there is only one sample of this sample type (i.e. sample->common_name)
  my $match = 0;
  my $sample;
  my $scn;
  foreach my $s (@dna_samples){
    my $current_scn = $s->common_name || "NULL";
    if ($tissue_type =~ /$current_scn/i){
      $match++;
      $sample = $s;
      $scn = $current_scn;
    }
  }
  if ($match == 0){
    $self->error_message("Did not find a matching DNA sample for tissue type: $tissue_type");
    exit 1;
  }elsif ($match > 1){
    $self->error_message("Found more than one matching DNA sample of tissue type: $tissue_type");
    exit 1;
  }else{
    $self->status_message("\nFound a DNA sample " . $sample->name . " ($scn) matching tissue type: $tissue_type");
  }

  #Is there WGS or Exome data?  Get the instrument data of each type for this sample.
  #- return if there is no data of the desired type
  my @tmp;
  my @sample_instrument_data = $self->get_dna_instrument_data('-sample'=>$sample, '-data_type'=>$data_type);
  return @tmp unless (scalar(@sample_instrument_data));

  my $subject_id = $sample->patient->id;
  my @models = $sample->models;
  my $model_count = scalar(@models);
  $self->status_message("\tStarting with " . $model_count . " models for this sample. Candidates that meet basic criteria:");

  my $ignore_models_bx;
  if (my $ignore_models_matching = $self->ignore_models_matching) {
    $ignore_models_bx = UR::BoolExpr->resolve_for_string('Genome::Model::ReferenceAlignment',$ignore_models_matching);
    $self->status_message("\tIgnoring models matching: $ignore_models_bx");
  }

  #Test for correct processing profile, reference sequence build, annotation build, and dbsnp build
  #Also make sure that all the instrument data of wgs or exome type is being used (exome can be exome+wgs lanes)
  foreach my $model (@models){
    next unless ($model->class eq "Genome::Model::ReferenceAlignment");
    next unless ($model->processing_profile_id eq $self->ref_align_pp->id);
    next unless ($model->reference_sequence_build->id eq $self->reference_sequence_build->id);
    if ($model->can("annotation_reference_build")){
      next unless ($model->annotation_reference_build);
      next unless ($model->annotation_reference_build->id eq $self->annotation_build->id);
    }else{
      next;
    }
    if ($model->can("dbsnp_build")){
      my $dbsnp_build = $model->dbsnp_build;
      if ($dbsnp_build){
        next unless ($dbsnp_build->id eq $self->dbsnp_build->id);
      }else{
        next;
      }
    }else{
      next;
    }
    if ($ignore_models_bx and $ignore_models_bx->evaluate($model)) {
        $self->status_message("\t\tIgnore: " . $model->__display_name__);
        next;
    }

    #Get genotype microarray model for this sample
    my $genotype_microarray_model_id = $self->get_genotype_microarray_model_id('-sample'=>$sample);

    #If genotype microarray data is available:
    #Make sure a genotype microarray model is defined for this ref-align model
    #If it is defined, make sure it is the correct one
    if ($sample->default_genotype_data){
      next unless ($model->genotype_microarray_model);
      next unless ($model->genotype_microarray_model->id eq $genotype_microarray_model_id);
    }

    #Is this a WGS or an Exome model?
    #WGS models do not have exome data, Exome models have at least one lane of exome data - skip those that are not the current type being considered
    next unless ($self->determine_model_data_type('-model'=>$model) eq $data_type);

    #If the desired $data_type is exome.  Check that the TRSN and ROI have been set correctly, exclude models that are not
    next unless ($self->check_model_trsn_and_roi('-model'=>$model, '-data_type'=>$data_type));
    #$self->status_message("\t\tName: " . $model->name . " (" . $model->id . ")");

    push (@final_models, $model);
  }

  #Make sure all the wgs or exome data is associated with the model
  #In both wgs and exome models, additional data will be allowed to handle weird situations.  
  @final_models = @{$self->check_for_missing_and_excluded_data('-models'=>\@final_models, '-sample_instrument_data'=>\@sample_instrument_data)};

  my $final_model_count = scalar(@final_models);
  $self->status_message("\tFound " . $final_model_count . " suitable $data_type models (matching default or user specified criteria):");
  foreach my $model (@final_models){
    $self->status_message("\t\tName: " . $model->name . " (" . $model->id . ")");
  }

  #If there are no suitable models, one will need to be created
  unless ($final_model_count > 0){
    $self->create_ref_align_model('-data_type'=>$data_type, '-sample'=>$sample, '-sample_instrument_data'=>\@sample_instrument_data);
    return @tmp;
  }

  #If there are suitable models, check the status of their *builds*, and if neccessary launch a new build
  my $models_status = $self->check_models_status('-models'=>\@final_models);

  #If there is one or more suitable successful models, return the models objects.  Must return all suitable models to allow checking against available somatic variation models
  if ($models_status){
    return @final_models;
  }else{
    return @tmp;
  }
}


#Gather rnaseq models for a single subject that meet all the specified criteria
sub check_rnaseq_models{
  my $self = shift;
  my %args = @_;
  my $tissue_type = $args{'-tissue_type'};
  my @rna_samples = @{$args{'-rna_samples'}};

  my @final_models;
  my @tmp;

  #Make sure there is only one sample of this sample type (i.e. sample->common_name)
  my $match = 0;
  my $sample;
  my $scn;
  foreach my $s (@rna_samples){
    my $current_scn = $s->common_name || "NULL";
    if ($current_scn =~ /$tissue_type/i){
      $match++;
      $sample = $s;
      $scn = $current_scn;
    }
  }
  if ($match == 0 && ($tissue_type =~ /normal/i)){
    #$self->warning_message("Did not find a matching RNA sample for tissue type: $tissue_type");
    return @tmp;
  }elsif($match == 0){
    $self->warning_message("Did not find a matching RNA sample for tissue type: $tissue_type");
    return @tmp;    
  }elsif ($match > 1){
    $self->error_message("Found more than one matching RNA sample of tissue type: $tissue_type");
  }else{
    $self->status_message("\nFound an RNA sample " . $sample->name . " ($scn) matching tissue type: $tissue_type");
  }

  #Is there actually any RNA-seq data?
  #- return if there is no data of the desired type
  my @sample_instrument_data = $sample->instrument_data;
  @sample_instrument_data = @{$self->exclude_instrument_data('-instrument_data'=>\@sample_instrument_data)};

  unless (scalar(@sample_instrument_data)){
    $self->status_message("\tCould not find any rna-seq data");
    return @tmp;
  }
  my $subject_id = $sample->patient->id;
  my @models = $sample->models;
  my $model_count = scalar(@models);
  $self->status_message("\tStarting with " . $model_count . " models for this sample. Candidates that meet basic criteria:");

  #Test for correct processing profile, reference sequence build, annotation build
  foreach my $model (@models){
    next unless ($model->class eq "Genome::Model::RnaSeq");
    next unless ($model->processing_profile_id eq $self->rnaseq_pp->id);
    next unless ($model->reference_sequence_build->id eq $self->reference_sequence_build->id);
    next unless ($model->annotation_build->id eq $self->annotation_build->id);
    push (@final_models, $model);
    #$self->status_message("\t\tName: " . $model->name . " (" . $model->id . ")");
  }

  #Make sure all the rna-seq data is associated with the model
  @final_models = @{$self->check_for_missing_and_excluded_data('-models'=>\@final_models, '-sample_instrument_data'=>\@sample_instrument_data)};

  my $final_model_count = scalar(@final_models);
  $self->status_message("\tFound " . $final_model_count . " suitable $tissue_type rna-seq models (matching default or user specified criteria):");
  foreach my $model (@final_models){
    $self->status_message("\t\tName: " . $model->name . " (" . $model->id . ")");
  }

  #If there are no suitable models, one will need to be created
  unless ($final_model_count > 0){
    $self->create_rnaseq_model('-sample'=>$sample, '-sample_instrument_data'=>\@sample_instrument_data);
    return @tmp;
  }

  #If there are suitable models, check the status of their *builds*, and if neccessary launch a new build
  my $models_status = $self->check_models_status('-models'=>\@final_models);

  #If there is one or more suitable successful models, return the models objects.  Must return all suitable models to allow checking against available somatic variation models
  if ($models_status){
    return @final_models;
  }else{
    return @tmp;
  }
}


#Check for existence of a differential expression model meeting desired criteria and create one if it does not exist
sub check_de_models{
  my $self = shift;
  my %args = @_;
  my $normal_rnaseq_model = $args{'-normal_rnaseq_model'};
  my $tumor_rnaseq_model = $args{'-tumor_rnaseq_model'};

  my @final_models;
  my @tmp;

  my @models = $self->differential_expression_pp->models;

  #Test for correct processing profile and input models
  foreach my $model (@models){
    next unless ($model->class eq "Genome::Model::DifferentialExpression");
    next unless ($model->processing_profile_id eq $self->differential_expression_pp->id);
    my $condition_model_ids_string = $model->condition_model_ids_string;
    my @groups = split(" ", $condition_model_ids_string);
    next unless (scalar(@groups) == 2);
    my @group1_members = split(",", $groups[0]);
    my @group2_members = split(",", $groups[1]);
    next unless (scalar(@group1_members == 1) && scalar(@group2_members));

    next unless ($group1_members[0] eq $normal_rnaseq_model->id);
    next unless ($group2_members[0] eq $tumor_rnaseq_model->id);

    #$self->status_message("\t\tName: " . $model->name . " (" . $model->id . ")");
    push (@final_models, $model);
  }

  my $final_model_count = scalar(@final_models);
  $self->status_message("\tFound " . $final_model_count . " suitable differential expression models (matching default or user specified criteria):");
  foreach my $model (@final_models){
    $self->status_message("\t\tName: " . $model->name . " (" . $model->id . ")");
  }

  #If there are no suitable models, one will need to be created
  unless ($final_model_count > 0){
    $self->create_de_model('-normal_rnaseq_model'=>$normal_rnaseq_model, '-tumor_rnaseq_model'=>$tumor_rnaseq_model);
    return @tmp;
  }

  my $models_status = $self->check_models_status('-models'=>\@final_models);

  #If there is one or more suitable successful models, return the models objects.  Must return all suitable models to allow checking against available clinseq models
  if ($models_status){
    return @final_models;
  }else{
    return @tmp;
  }
}


#Check for existence of a somatic variation model meeting desired criteria and create one if it does not exist
sub check_somatic_variation_models{
  my $self = shift;
  my %args = @_;
  my $data_type = $args{'-data_type'};
  my @dna_samples = @{$args{'-dna_samples'}};
  my @normal_ref_align_models = @{$args{'-normal_models'}};
  my @tumor_ref_align_models = @{$args{'-tumor_models'}};

  my @final_models;
  my @tmp;

  my $somatic_variation_pp_id;
  if ($data_type eq 'wgs'){
    $somatic_variation_pp_id = $self->wgs_somatic_variation_pp->id;
  }elsif($data_type eq 'exome'){
    $somatic_variation_pp_id = $self->exome_somatic_variation_pp->id;
  }else{
    $self->error_message("Data type must be wgs or exome but $data_type was provided");
    exit 1;
  }

  #Get models for either tumor or normal DNA sample 
  my @existing_models;
  foreach my $sample(@dna_samples){
    my @models = $sample->models;
    push(@existing_models, @models);
  }
  my $existing_model_count = scalar(@existing_models);
  $self->status_message("\tStarting with " . $existing_model_count . " models for this pair of DNA samples. Candidates that meet basic criteria:");

  my $ignore_models_bx;
  if (my $ignore_models_matching = $self->ignore_models_matching) {
    $ignore_models_bx = UR::BoolExpr->resolve_for_string('Genome::Model::SomaticVariation',$ignore_models_matching);
    $self->status_message("\tIgnoring models matching: $ignore_models_bx");
  }

  #Test for correct processing profile, annotation build
  foreach my $model (@existing_models){
    next unless ($model->class eq "Genome::Model::SomaticVariation");
    next unless ($model->processing_profile_id eq $somatic_variation_pp_id);
    next unless ($model->annotation_build->id eq $self->annotation_build->id);
    if ($ignore_models_bx and $ignore_models_bx->evaluate($model)) {
        $self->status_message("\t\tIgnore: " . $model->__display_name__);
        next;
    }

    #If previously discovered variants was defined as an input to the somatic-variation model, disallow it...
    #-> Don't do this now that previously discovered variations build is only used to annotate variants with dbSNP ids
    #if ($model->can("previously_discovered_variations_build")){
    #  next if ($model->previously_discovered_variations_build);
    #}

    #Check for the correct version of previously discovered variants
    my $got = $model->previously_discovered_variations_build;
    my $expected = $self->previously_discovered_variations;
    unless ($got){
      $self->status_message("skipping " . $model->__display_name__ . " because it does not have previously discovered variants set");
      next;
    }
    unless ($got->id eq $expected->id){
      $self->status_message("\tskipping " . $model->__display_name__ . " because previously discovered variations build "
        . $got->__display_name__ 
        . " doesn't match expected " . $expected->__display_name__
      );
      next;
    }
    
    #Make sure one of the passing normal AND tumor reference alignment models are specified as inputs to the somatic variation model
    my $tumor_model_id = $model->tumor_model->id;
    my $tumor_model_match = 0;
    foreach my $model (@tumor_ref_align_models){
      $tumor_model_match = 1 if ($model->id eq $tumor_model_id);
    }
    my $normal_model_id = $model->normal_model->id;
    my $normal_model_match = 0;
    foreach my $model (@normal_ref_align_models){
      $normal_model_match = 1 if ($model->id eq $normal_model_id);
    }
    next unless ($tumor_model_match && $normal_model_match);

    #$self->status_message("\t\tUse: " . $model->__display_name__);
    push (@final_models, $model);
  }

  my $final_model_count = scalar(@final_models);
  $self->status_message("\tFound " . $final_model_count . " suitable $data_type somatic-variation models (matching default or user specified criteria):");
  foreach my $model (@final_models){
    $self->status_message("\t\tName: " . $model->name . " (" . $model->id . ")");
  }

  #If there are no suitable models, one will need to be created
  unless ($final_model_count > 0){
    $self->create_somatic_variation_model('-processing_profile_id'=>$somatic_variation_pp_id, '-normal_models'=>\@normal_ref_align_models, '-tumor_models'=>\@tumor_ref_align_models);
    return @tmp;
  }

  #TODO: Now check the input builds of supposedly good somatic-variation models to make sure they are using the latest reference-alignment builds
  #We already check to see if ref-align data needs to be added to a ref-align model...  but if that does happen the somatic-models are out-of-date
  @final_models = @{$self->check_somatic_input_builds('-models'=>\@final_models)};
  $final_model_count = scalar(@final_models);
  return @tmp unless ($final_model_count > 0);

  #If there are suitable models, check the status of their *builds*, and if neccessary launch a new build
  my $models_status = $self->check_models_status('-models'=>\@final_models);

  #If there is one or more suitable successful models, return the models objects.  Must return all suitable models to allow checking against available clinseq models
  if ($models_status){
    return @final_models;
  }else{
    return @tmp;
  }

  return;
}

sub check_somatic_input_builds{
  my $self = shift;
  my %args = @_;
  my @models = @{$args{'-models'}};
  my @final_models;

  foreach my $somatic_model (@models){
    my $latest_somatic_build = $somatic_model->latest_build;
    my $last_complete_somatic_build = $somatic_model->last_succeeded_build;

    my $latest_tumor_build = $latest_somatic_build->tumor_build;
    my $latest_normal_build = $latest_somatic_build->normal_build;

    my $tumor_model = $somatic_model->tumor_model;
    my $lc_tumor_build = $tumor_model->last_succeeded_build;
    my $normal_model = $somatic_model->normal_model;
    my $lc_normal_build = $normal_model->last_succeeded_build;

    my $somatic_model_id = $somatic_model->id;

    if ($last_complete_somatic_build && $lc_tumor_build && $lc_normal_build){
      #There is at least one complete somatic build
      my $last_complete_tumor_build = $last_complete_somatic_build->tumor_build;
      my $last_complete_normal_build = $last_complete_somatic_build->normal_build;
      if (($last_complete_tumor_build->id eq $lc_tumor_build->id) && ($last_complete_normal_build->id eq $lc_normal_build->id)){
        push(@final_models, $somatic_model);
      }elsif(($latest_tumor_build->id eq $lc_tumor_build->id) && ($latest_normal_build->id eq $lc_normal_build->id)){
        $self->status_message("WARNING -> latest build of somatic model $somatic_model_id is using the latest refalign builds but has the following status: " . $latest_somatic_build->status);
      }else{
        $self->status_message("WARNING -> latest build of somatic model: $somatic_model_id is not using the last complete builds of the underlying refalign models.  Run the following:");
        $self->status_message("genome model build start $somatic_model_id");
      }
    }elsif($latest_somatic_build && $lc_tumor_build && $lc_normal_build){
      if(($latest_tumor_build->id eq $lc_tumor_build->id) && ($latest_normal_build->id eq $lc_normal_build->id)){
        $self->status_message("WARNING -> latest build of somatic model $somatic_model_id is using the latest refalign builds but has the following status: " . $latest_somatic_build->status);
      }else{
        $self->status_message("WARNING -> latest build of somatic model: $somatic_model_id is not using the last complete builds of the underlying refalign models.  Run the following:");
        $self->status_message("genome model build start $somatic_model_id");
      }
    }
  }

  return(\@final_models);
}

#Create a reference-alignment model for wgs or exome data for a single sample
sub create_ref_align_model{
  my $self = shift;
  my %args = @_;
  my $data_type = $args{'-data_type'};
  my $sample = $args{'-sample'};
  my @sample_instrument_data = @{$args{'-sample_instrument_data'}}; #Already limited to data of $data_type
 
  #Use the predefined list of instrument data instead of the '--all' option.
  #This will help if we want to exclude instrument data for some models
  my @iids;
  foreach my $instrument_data (@sample_instrument_data){
    push (@iids, $instrument_data->id);
  }
  my $iids_list = join(",", @iids);

  #Get genotype microarray model for this sample
  my $genotype_microarray_model_id = $self->get_genotype_microarray_model_id('-sample'=>$sample);

  #If genotype microarray data is available but a successful model is not, do not proceed with creation of reference alignment models
  if ($sample->default_genotype_data && $genotype_microarray_model_id eq "0"){
    $self->status_message("WARNING -> genotype microarray data was found but a model with the desired criteria could not be found.  Check for issues with your existing microarray models");
    return;
  }

  #Come up with a descriptive model name.  May be more trouble than it is worth, perhaps it should be autogenerated?
  my $sample_name = $sample->name;
  my $annotation_id = $self->annotation_build->id;
  my $reference_build_id = $self->reference_sequence_build->id;
  my $ref_align_pp_id = $self->ref_align_pp->id;
  my $dbsnp_build_id = $self->dbsnp_build->id;

  #Only specify a genotype microarray model id as input, if one was found
  my $genotype_microarray_string = '';
  if ($genotype_microarray_model_id){
    $genotype_microarray_string = "--genotype-microarray-model='$genotype_microarray_model_id'";
  }

  my @commands;

  #WGS data
  if ($data_type eq 'wgs'){
    push(@commands, "\n#Create a WGS reference-alignment model as follows:");
    push(@commands, "genome model define reference-alignment  --reference-sequence-build='$reference_build_id'  --annotation-reference-build='$annotation_id'  --subject='$sample_name'  --processing-profile='$ref_align_pp_id'  --dbsnp-build='$dbsnp_build_id'  --instrument-data='$iids_list'  $genotype_microarray_string");
    push(@commands, "genome model build start ''");
  }

  #Exome data
  if ($data_type eq 'exome'){
    #Determine the target region set name associated with this list of instrument-data
    my $target_region_set_name = $self->get_trsn('-instrument_data'=>\@sample_instrument_data);

    #Determine the desired region of interest set name for this target region set name
    my $region_of_interest_set_name = $self->get_roi_name('-target_region_set_name'=>$target_region_set_name);

    push(@commands, "\n#Create an Exome reference-alignment model as follows:");
    push(@commands, "genome model define reference-alignment  --reference-sequence-build='$reference_build_id'  --annotation-reference-build='$annotation_id'  --subject='$sample_name'  --processing-profile='$ref_align_pp_id'  --dbsnp-build='$dbsnp_build_id'  --target-region-set-names='$target_region_set_name'  --region-of-interest-set-name='$region_of_interest_set_name'  --instrument-data='$iids_list'  $genotype_microarray_string");
    push(@commands, "genome model build start ''");
  }
  foreach my $line (@commands){
    $self->status_message($line);
  }

  return;
}


#Create a rna-seq model for rnaseq data for a single sample
sub create_rnaseq_model{
  my $self = shift;
  my %args = @_;
  my $sample = $args{'-sample'};
  my @sample_instrument_data = @{$args{'-sample_instrument_data'}}; #Already limited to data of $data_type
 
  #Use the predefined list of instrument data instead of the '--all' option.
  #This will help if we want to exclude instrument data for some models
  my @iids;
  foreach my $instrument_data (@sample_instrument_data){
    push (@iids, $instrument_data->id);
  }
  my $iids_list = join(",", @iids);

  my $sample_name = $sample->name;
  my $annotation_id = $self->annotation_build->id;
  my $reference_build_id = $self->reference_sequence_build->id;
  my $rnaseq_pp_id = $self->rnaseq_pp->id;

  my @commands;

  push(@commands, "\n#Create an RNA-seq model as follows:");
  push(@commands, "genome model define rna-seq  --reference-sequence-build='$reference_build_id'  --annotation-build='$annotation_id'  --subject='$sample_name'  --processing-profile='$rnaseq_pp_id'  --instrument-data='$iids_list'");
  push(@commands, "genome model build start ''");

  foreach my $line (@commands){
    $self->status_message($line);
  }

  return;
}


#Create a differential-expression model for a singel tumor/normal rna-seq pair
sub create_de_model{
  my $self = shift;
  my %args = @_;
  my $normal_rnaseq_model = $args{'-normal_rnaseq_model'};
  my $tumor_rnaseq_model = $args{'-tumor_rnaseq_model'};

  my $individual_id = $self->individual->id;
  my $individual_name = $self->individual->name;
  my $normal_rnaseq_model_id = $normal_rnaseq_model->id;
  my $normal_subject_name = $normal_rnaseq_model->subject->name;
  my $normal_individual_common_name = $normal_rnaseq_model->subject->patient_common_name;
  my $normal_sample_common_name = $normal_rnaseq_model->subject->common_name;
  my $tumor_rnaseq_model_id = $tumor_rnaseq_model->id;
  my $tumor_subject_name = $tumor_rnaseq_model->subject->name;
  my $tumor_individual_common_name = $tumor_rnaseq_model->subject->patient_common_name;
  my $tumor_sample_common_name = $tumor_rnaseq_model->subject->common_name;
  my $differential_expression_pp_id = $self->differential_expression_pp->id;

  my $final_normal_name = $normal_subject_name;
  my $final_tumor_name = $tumor_subject_name;

  if ($normal_sample_common_name && $tumor_sample_common_name){
    unless ($normal_sample_common_name eq $tumor_sample_common_name){
      $final_normal_name = $normal_sample_common_name;
      $final_tumor_name = $tumor_sample_common_name;
    }
  }

  #Remove trailing and leading spaces
  $final_normal_name =~ s/^\s+//; 
  $final_normal_name =~ s/\s+$//;
  $final_tumor_name =~ s/^\s+//; 
  $final_tumor_name =~ s/\s+$//;

  my @commands;

  #genome model define differential-expression --processing-profile=? --condition-model-ids-string=? --condition-labels-string=? 
  push(@commands, "\n#Create a differential-expression model as follows:");
  push(@commands, "genome model define differential-expression --processing-profile='$differential_expression_pp_id' --condition-model-ids-string='$normal_rnaseq_model_id $tumor_rnaseq_model_id' --condition-labels-string='$final_normal_name,$final_tumor_name' --subject='$individual_name'");
  push(@commands, "genome model build start ''");

  foreach my $line (@commands){
    $self->status_message($line);
  }

  return;
}


#Create a somatic-variation model for a single tumor/normal pair
sub create_somatic_variation_model{
  my $self = shift;
  my %args = @_;
  my $processing_profile_id = $args{'-processing_profile_id'};
  my @normal_models = @{$args{'-normal_models'}};
  my @tumor_models = @{$args{'-tumor_models'}};

  #If there are multiple suitable normal or tumor reference alignment models, select the 'best' one
  #The definition of 'best' is not neccessarily clear... Try selecting that with the most data, break ties by creation data (favoring recently created models)
  my $best_tumor_model = $self->select_best_model('-models'=>\@tumor_models);
  my $best_normal_model = $self->select_best_model('-models'=>\@normal_models);

  my $tumor_model_id = $best_tumor_model->id;
  my $normal_model_id = $best_normal_model->id;
  my $somatic_variation_pp_id = $processing_profile_id;
  my $annotation_build_id = $self->annotation_build->id;
  my $previously_discovered_variations_build_id = $self->previously_discovered_variations->id;

  #Make sure neither of the input models/builds has been archived before proceeding...
  my $tumor_build = $best_tumor_model->last_succeeded_build;
  unless ($self->skip_check_archived){
    if ($tumor_build->is_archived){
      $self->status_message("\tWARNING -> Tumor build is currently archived. Run: bsub genome model build unarchive --lab=Mardis-Wilson " . $tumor_build->id);
      return;
    }
  }
  my $normal_build = $best_normal_model->last_succeeded_build;
  unless ($self->skip_check_archived){
    if ($normal_build->is_archived){
      $self->status_message("\tWARNING -> Normal build is currently archived. Run: bsub genome model build unarchive --lab=Mardis-Wilson " . $normal_build->id);
      return;
    }
  }
  my @commands;
  push(@commands, "\n#Create a Somatic-Variation model as follows:");
  push(@commands, "genome model define somatic-variation  --processing-profile=$processing_profile_id  --tumor-model=$tumor_model_id  --normal-model=$normal_model_id  --annotation-build=$annotation_build_id  --previously-discovered-variations-build=$previously_discovered_variations_build_id");
  push(@commands, "genome model build start ''");

  foreach my $line (@commands){
    $self->status_message($line);
  }

  return;
}


#Select the best model from a list of suitable models according to simple criteria that may depend on the model type being considered
sub select_best_model{
  my $self = shift;
  my %args = @_;
  my @models = @{$args{'-models'}};

  #If there are no models, return 0
  return 0 if (scalar(@models) == 0);

  #If there is only one model, it is the best by default
  my $best_model = $models[0];
  return $best_model if (scalar(@models) == 1);

  #Gather ranking information for the models
  my %candidate_models;
  my $m = 0;
  my $max_id_count = 0;
  foreach my $model (@models){
    $m++;
    $candidate_models{$m}{model} = $model;
    my $id_count = 0;
    if ($model->can("instrument_data")){
      my @instrument_data = $model->instrument_data;
      $id_count = scalar(@instrument_data);
    }
    $candidate_models{$m}{id_count} = $id_count;
    $max_id_count = $id_count if ($id_count > $max_id_count);
  }

  #Rank models according to amount of instrument data available
  foreach my $m (keys %candidate_models){
    $best_model = $candidate_models{$m}{model};
    delete $candidate_models{$m} unless ($candidate_models{$m}{id_count} == $max_id_count);
  }
  my $candidate_model_count = keys %candidate_models;
  return $best_model if ($candidate_model_count == 1);

  #If there is still no clear winner, rank according to highest id of model
  my $max_model_id = 0;
  my $latest_date = "1900-01-01 00:00:00";
  my $dateformat = "%Y-%m-%d %H:%M:%S";

  foreach my $m (keys %candidate_models){
    my $model = $candidate_models{$m}{model};
    my $model_creation_date = $model->creation_date;
    my $date1 = Time::Piece->strptime($latest_date, $dateformat);
    my $date2 = Time::Piece->strptime($model_creation_date, $dateformat);
    if ($date2 > $date1){
      $latest_date = $model_creation_date;
      $best_model = $model;
    }
  }
  return $best_model;
}


#Check model build status for each model.  Return true if it is safe to proceed, launch builds where needed, report if we need to wait for running builds
sub check_models_status{
  my $self = shift;
  my %args = @_;
  my @models = @{$args{'-models'}};
  my $silent = 0;
  if ($args{'-silent'}){
    $silent = $args{'-silent'};
  }

  #Foreach model:
  #- if there is a successful build do nothing for this model
  #- if there is a running or scheduled build, do nothing for this model
  #- if there is an unstartable build warn the user?
  #- if there are no builds at all or only failed builds and abandoned builds (i.e. if the previous scenarios are not true), start one 
  my $model_count = scalar(@models);
  my $ready_model_count = 0;
  foreach my $model (@models){
    my $model_id = $model->id;
    my @builds = $model->builds;
    my $build_count = scalar(@models);
  
    my $running_builds = 0;
    my $succeeded_builds = 0;
    my $unstartable_builds = 0;
    my $scheduled_builds = 0;
    foreach my $build (@builds){
      my $build_id = $build->id;
      my $status = $build->status;
      if ($status =~ /succeeded/i){
        $succeeded_builds++;
      }elsif ($status =~ /running/i){
        $running_builds++;
      }elsif ($status =~ /unstartable/i){
        $unstartable_builds++;
      }elsif ($status =~ /scheduled/i){
        $scheduled_builds++;
      }
    }

    #For each model, get the last succeeded build and check if it could be set as do-not-archive
    #Before attempting to use a model/build as input to another model, make sure it is not archived (and check whether it has been set as do-not-archive)
    #- Do not start automatically marking these things as do-not-archive, leave it to the end user to mark the actual desired clin-seq models as do-not-archive
    foreach my $model (@models){
      my $model_id = $model->id;
      my $build = $model->last_succeeded_build;
      if ($build){
        my $build_id = $build->id;
        unless ($self->skip_check_archived){
          if ($build->is_archived){
            $self->status_message("\tWARNING -> Successful build $build_id of model $model_id that meets desired criteria is currently archived! Consider running: bsub genome model build unarchive --lab=Mardis-Wilson $build_id") unless $silent;
          }elsif ($build->archivable){
            if ($self->check_archivable_status){
              $self->status_message("\tSuccessful build $build_id of model $model_id that meets desired criteria is currently archivable. Consider running: genome model build set-do-not-archive --reason='build needed for clin-seq model' $build_id") unless $silent;
            }
          }
        }
      }
    }

    #If a model has unstartable builds and now running or succeeded builds, warn the user
    if ($unstartable_builds && $running_builds == 0 && $succeeded_builds == 0){
      $self->status_message("\tWARNING -> Model: $model_id has unstartable builds only! Please investigate...") unless $silent;
    }
    
    if ($running_builds){
      $self->status_message("\tWARNING -> Model: $model_id has running builds...") unless $silent;
    }elsif ($succeeded_builds){
      $ready_model_count++;
    }elsif ($scheduled_builds){
      $self->status_message("\tWARNING -> Model: $model_id has a scheduled build...") unless $silent;
    }else{
      $self->status_message("\tWARNING -> Model: $model_id needs a build ...\n\tgenome model build start $model_id") unless $silent;
    }
  }

  #If all models are 'ready' return true, or perhaps if at least one model is 'ready' return true
  if ($ready_model_count > 0 && $ready_model_count == $model_count){
    return 1;
  }
  return undef;
}


#Check status of existing clinseq models with a specified list of input models and create one if needed
sub check_clinseq_models{
  my $self = shift;
  my %args = @_;
  my $wgs_model = $args{'-wgs_somatic_variation_model'};
  my $exome_model = $args{'-exome_somatic_variation_model'};
  my $normal_rnaseq_model = $args{'-normal_rnaseq_model'};
  my $tumor_rnaseq_model = $args{'-tumor_rnaseq_model'};
  my $de_model = $args{'-de_model'};

  my @samples = @{$args{'-samples'}};

  my @final_models;
  my @tmp;

  #If there are no candidate input models at all, skip
  return @tmp unless ($wgs_model || $exome_model || $normal_rnaseq_model || $tumor_rnaseq_model);

  #Based on the types of data available for this patient do we have all the model types we hope to have?
  #i.e. if there is all four data types available for the sample, there should be all four input models
  #Warn if not, but attempt to check or create a more limited clin-seq model anyway... (e.g. wgs & exome somatic variation but not rna-seq)
  my $wgs_data = $self->check_instrument_data('-data_type'=>'wgs', '-samples'=>\@samples);
  my $exome_data = $self->check_instrument_data('-data_type'=>'exome', '-samples'=>\@samples);
  my $normal_rnaseq_data = $self->check_instrument_data('-data_type'=>'normal_rnaseq', '-samples'=>\@samples);
  my $tumor_rnaseq_data = $self->check_instrument_data('-data_type'=>'tumor_rnaseq', '-samples'=>\@samples);

  if ($wgs_data == 2 && !$wgs_model){
    $self->warning_message("WGS data exists but no suitable, complete WGS somatic-variation model was found")
  }
  if ($exome_data == 2 && !$exome_model){
    $self->warning_message("Exome data exists but no suitable, complete Exome somatic-variation model was found");
  }
  if ($normal_rnaseq_data && !$normal_rnaseq_model){
    $self->warning_message("Normal RNAseq data exists but no suitable, complete normal RNAseq model was found");
  }
  if ($tumor_rnaseq_data && !$tumor_rnaseq_model){
    $self->warning_message("Tumor RNAseq data exists but no suitable, complete tumor RNAseq model was found");
  }
  if (($normal_rnaseq_model && $tumor_rnaseq_model) && !$de_model){
    $self->warning_message("Normal and Tumor RNAseq models exist, but no suitable, complete differential-expression model was found");
  }

  my @existing_models = Genome::Model->get(processing_profile_id => $self->clinseq_pp->id, subject_id => $self->individual->id);
  my $existing_model_count = scalar(@existing_models);
  $self->status_message("\tStarting with " . $existing_model_count . " clin-seq models for this individual. Candidates that meet basic criteria:");

  #Create a hash of the samples specified by the user.  
  my %target_samples;
  foreach my $sample (@samples){
    $target_samples{$sample->id}=1;
  }

  #Make sure the 'best' input models of each type are being used
  foreach my $model (@existing_models){
    my $current_wgs_model = $model->wgs_model;
    my $current_exome_model = $model->exome_model;
    my $current_normal_rnaseq_model = $model->normal_rnaseq_model;
    my $current_tumor_rnaseq_model = $model->tumor_rnaseq_model;
    my $current_de_model = $model->de_model;

    #Only consider clin-seq models whose input models are using the correct processing profiles
    if ($current_wgs_model){
      next unless ($current_wgs_model->processing_profile->id eq $self->wgs_somatic_variation_pp->id);
    }
    if ($current_exome_model){
      next unless ($current_exome_model->processing_profile->id eq $self->exome_somatic_variation_pp->id);
    }
    if ($current_normal_rnaseq_model){
      next unless ($current_normal_rnaseq_model->processing_profile->id eq $self->rnaseq_pp->id);
    }
    if ($current_tumor_rnaseq_model){
      next unless ($current_tumor_rnaseq_model->processing_profile->id eq $self->rnaseq_pp->id);
    }
    if ($current_de_model){
      next unless ($current_de_model->processing_profile->id eq $self->differential_expression_pp->id);
    }

    #If any sample on any of the input models is not in the user supplied list of target samples, that clin-seq model is not valid
    my %model_samples;
    my $found_nonmatching_sample = 0;
    if ($current_wgs_model){
      $model_samples{$current_wgs_model->normal_model->subject->id}=1;
      $model_samples{$current_wgs_model->tumor_model->subject->id}=1;
    }
    if ($current_exome_model){
      $model_samples{$current_exome_model->normal_model->subject->id}=1;
      $model_samples{$current_exome_model->tumor_model->subject->id}=1;
    }
    if ($current_normal_rnaseq_model){
      $model_samples{$current_normal_rnaseq_model->subject->id}=1;
    }
    if ($current_tumor_rnaseq_model){
      $model_samples{$current_tumor_rnaseq_model->subject->id}=1;
    }
    foreach my $current_sample_id (keys %model_samples){
      unless ($target_samples{$current_sample_id}){
        $found_nonmatching_sample = 1;
      }
    }
    next if ($found_nonmatching_sample);

    #Only consider clin-seq models whose annotation inputs match the current defaults defined for the pipeline
    my $cancer_annotation_db_id = Genome::Model::ClinSeq->__meta__->property("cancer_annotation_db")->default_value;
    my $cancer_annotation_db = Genome::Db::Tgi::CancerAnnotation->get($cancer_annotation_db_id);
    next unless ($model->cancer_annotation_db->id eq $cancer_annotation_db_id);

    my $misc_annotation_db_id = Genome::Model::ClinSeq->__meta__->property("misc_annotation_db")->default_value;
    my $misc_annotation_db = Genome::Db::Tgi::MiscAnnotation->get($misc_annotation_db_id);
    next unless ($model->misc_annotation_db->id eq $misc_annotation_db_id);

    my $cosmic_annotation_db_id = Genome::Model::ClinSeq->__meta__->property("cosmic_annotation_db")->default_value;
    my $cosmic_annotation_db = Genome::Db::Cosmic->get($cosmic_annotation_db_id);
    next unless ($model->cosmic_annotation_db->id eq $cosmic_annotation_db_id);

    #If a 'best' model is defined, only compare against clinseq models that have a model of that type defined and skip if the actual model does not match
    if ($wgs_model){
      next unless $current_wgs_model;
      next unless ($wgs_model->id eq $current_wgs_model->id);
    }
    if ($exome_model){
      next unless $current_exome_model;
      next unless ($exome_model->id eq $current_exome_model->id);
    }
    if ($normal_rnaseq_model){
      next unless $current_normal_rnaseq_model;
      next unless ($normal_rnaseq_model->id eq $current_normal_rnaseq_model->id);
    }
    if ($tumor_rnaseq_model){
      next unless $current_tumor_rnaseq_model;
      next unless ($tumor_rnaseq_model->id eq $current_tumor_rnaseq_model->id);
    }
    if ($de_model){
      next unless $current_de_model;
      next unless ($de_model->id eq $current_de_model->id);
    }
    push(@final_models, $model);
  }

  my $final_model_count = scalar(@final_models);
  $self->status_message("\tFound " . $final_model_count . " suitable clin-seq models (matching default or user specified criteria)");

  #If there are no suitable models, one will need to be created
  unless ($final_model_count > 0){
    $self->create_clinseq_model('-wgs_somatic_variation_model'=>$wgs_model,
                                '-exome_somatic_variation_model'=>$exome_model,
                                '-normal_rnaseq_model'=>$normal_rnaseq_model,
                                '-tumor_rnaseq_model'=>$tumor_rnaseq_model,
                                '-de_model'=>$de_model);
    return @tmp;
  }

  #If there are suitable models, check the status of their *builds*, and if neccessary launch a new build
  #Narrow the field to successful models to return as well
  my @final_models_succeeded2;
  foreach my $model (@final_models){
    my @test;
    push (@test, $model);
    my $models_status = $self->check_models_status('-models'=>\@test);
    push (@final_models_succeeded2, $model) if $models_status;
  }

  if (scalar(@final_models_succeeded2)){
    #Select the best clinseq model and return it
    my $best_clinseq_model = $self->select_best_model('-models'=>\@final_models_succeeded2);
    return $best_clinseq_model;
  }else{
    return 0;
  }
}


#Create a clin-seq model based on the available input models
sub create_clinseq_model{
  my $self = shift;
  my %args = @_;
  my $wgs_model = $args{'-wgs_somatic_variation_model'};
  my $exome_model = $args{'-exome_somatic_variation_model'};
  my $normal_rnaseq_model = $args{'-normal_rnaseq_model'};
  my $tumor_rnaseq_model = $args{'-tumor_rnaseq_model'};
  my $de_model = $args{'-de_model'};
  my $final_individual_name = $self->get_final_individual_name;

  my $clinseq_pp_id = $self->clinseq_pp->id;
  my $wgs_model_id = $wgs_model->id if $wgs_model;
  my $exome_model_id = $exome_model->id if $exome_model;
  my $normal_rnaseq_model_id = $normal_rnaseq_model->id if $normal_rnaseq_model;
  my $tumor_rnaseq_model_id = $tumor_rnaseq_model->id if $tumor_rnaseq_model;
  my $de_model_id = $de_model->id if $de_model;

  #Make sure a last_succeeded_build is defined for all input models
  #Make sure none of the input models/builds have been archived before proceeding...
  my @test_models = ($wgs_model, $exome_model, $normal_rnaseq_model, $tumor_rnaseq_model, $de_model);
  my $ready = 1;
  foreach my $model (@test_models){
    next unless $model;
    my $build = $model->last_succeeded_build;
    if ($build){
      unless ($self->skip_check_archived){
        if ($build->is_archived){
          $ready = 0;
          $self->status_message("\tWARNING -> Build is currently archived for model: " . $model->name . "\n\tRun: bsub genome model build unarchive --lab=Mardis-Wilson " . $build->id);
        }
      }
    }else{
      $ready = 0;
      $self->status_message("\tWARNING -> No successful build for clinseq input model: " . $model->name);
    }
  }

  #Unless all defined models have a successful build and are not archived, do not proceed
  unless ($ready){
    $self->status_message("\tWARNING -> Not ready to proceed with creation of clin-seq model, see warnings above");
    return;
  }

  #genome model define clin-seq  --processing-profile='November 2011 Clinical Sequencing'  --wgs-model='2882504846'  --exome-model='2882505032'  --tumor-rnaseq-model='2880794613'
  my @commands;
  my $clinseq_cmd = "genome model define clin-seq  --processing-profile='$clinseq_pp_id'";
  $clinseq_cmd .= "  --wgs-model='$wgs_model_id'" if $wgs_model;
  $clinseq_cmd .= "  --exome-model='$exome_model_id'" if $exome_model;
  $clinseq_cmd .= "  --normal-rnaseq-model='$normal_rnaseq_model_id'" if $normal_rnaseq_model;
  $clinseq_cmd .= "  --tumor-rnaseq-model='$tumor_rnaseq_model_id'" if $tumor_rnaseq_model;
  $clinseq_cmd .= "  --de-model='$de_model_id'" if $de_model;

  push(@commands, "\n#Create a Clin-Seq model for $final_individual_name as follows:");
  push(@commands, $clinseq_cmd);
  push(@commands, "genome model build start ''");

  foreach my $line (@commands){
    $self->status_message($line);
  }

  return;
}


sub check_clinseq_inputs{
  my $self = shift;
  my %args = @_;
  my $clinseq_model = $args{'-model'};

  my $clinseq_inputs_ok = 1;

  my $clinseq_model_id = $clinseq_model->id;
  my $clinseq_build = $clinseq_model->last_complete_build;
  if ($clinseq_build){

    #Check wgs input build used in clinseq build if defined
    if ($clinseq_model->wgs_model){
      my $wgs_model = $clinseq_model->wgs_model;
      my $wgs_build = $clinseq_build->wgs_build;
      my $lc_wgs_build = $wgs_model->last_complete_build;
      if ($lc_wgs_build){
        unless ($lc_wgs_build->id eq $wgs_build->id){
          $self->status_message("WARNING: last_complete_build of wgs model " . $wgs_model->id . " associated with clinseq model $clinseq_model_id is not being used.  Run the following:");
          $self->status_message("genome model build start $clinseq_model_id");
          $clinseq_inputs_ok = 0;
        }
      }else{
        $self->status_message("WARNING: wgs model " . $wgs_model->id . " associated with clinseq model $clinseq_model_id does not have a last_complete_build");
        $clinseq_inputs_ok = 0;
      }
    }

    #Check exome input build used in clinseq build if defined
    if ($clinseq_model->exome_model){
      my $exome_model = $clinseq_model->exome_model;
      my $exome_build = $clinseq_build->exome_build;
      my $lc_exome_build = $exome_model->last_complete_build;
      if ($lc_exome_build){
        unless ($lc_exome_build->id eq $exome_build->id){
          $self->status_message("WARNING: last_complete_build of exome model " . $exome_model->id . " associated with clinseq model $clinseq_model_id is not being used.  Run the following:");
          $self->status_message("genome model build start $clinseq_model_id");
          $clinseq_inputs_ok = 0;
        }
      }else{
        $self->status_message("WARNING: exome model " . $exome_model->id . " associated with clinseq model $clinseq_model_id does not have a last_complete_build");
        $clinseq_inputs_ok = 0;
      }
    }

    #Check normal rnaseq input build used in clinseq build if defined
    if ($clinseq_model->normal_rnaseq_model){
      my $normal_rnaseq_model = $clinseq_model->normal_rnaseq_model;
      my $normal_rnaseq_build = $clinseq_build->normal_rnaseq_build;
      my $lc_normal_rnaseq_build = $normal_rnaseq_model->last_complete_build;
      if ($lc_normal_rnaseq_build){
        unless ($lc_normal_rnaseq_build->id eq $normal_rnaseq_build->id){
          $self->status_message("WARNING: last_complete_build of normal_rnaseq model " . $normal_rnaseq_model->id . " associated with clinseq model $clinseq_model_id is not being used.  Run the following:");
          $self->status_message("genome model build start $clinseq_model_id");
          $clinseq_inputs_ok = 0;
        }
      }else{
        $self->status_message("WARNING: normal_rnaseq model " . $normal_rnaseq_model->id . " associated with clinseq model $clinseq_model_id does not have a last_complete_build");
        $clinseq_inputs_ok = 0;
      }
    }
    
    #Check tumor rnaseq input build used in clinseq build if defined
    if ($clinseq_model->tumor_rnaseq_model){
      my $tumor_rnaseq_model = $clinseq_model->tumor_rnaseq_model;
      my $tumor_rnaseq_build = $clinseq_build->tumor_rnaseq_build;
      my $lc_tumor_rnaseq_build = $tumor_rnaseq_model->last_complete_build;
      if ($lc_tumor_rnaseq_build){
        unless ($lc_tumor_rnaseq_build->id eq $tumor_rnaseq_build->id){
          $self->status_message("WARNING: last_complete_build of tumor_rnaseq model " . $tumor_rnaseq_model->id . " associated with clinseq model $clinseq_model_id is not being used.  Run the following:");
          $self->status_message("genome model build start $clinseq_model_id");
          $clinseq_inputs_ok = 0;
        }
      }else{
        $self->status_message("WARNING: tumor_rnaseq model " . $tumor_rnaseq_model->id . " associated with clinseq model $clinseq_model_id does not have a last_complete_build");
        $clinseq_inputs_ok = 0;
      }
    }

    #Check de input build used in clinseq build if defined
    if ($clinseq_model->de_model){
      my $de_model = $clinseq_model->de_model;
      my $de_build = $clinseq_build->de_build;
      my $lc_de_build = $de_model->last_complete_build;
      if ($lc_de_build){
        unless ($lc_de_build->id eq $de_build->id){
          $self->status_message("WARNING: last_complete_build of de model " . $de_model->id . " associated with clinseq model $clinseq_model_id is not being used.  Run the following:");
          $self->status_message("genome model build start $clinseq_model_id");
          $clinseq_inputs_ok = 0;
        }
      }else{
        $self->status_message("WARNING: de model " . $de_model->id . " associated with clinseq model $clinseq_model_id does not have a last_complete_build");
        $clinseq_inputs_ok = 0;
      }
    }

  }else{
    $self->status_message("WARNING: There is no last_complete_build for this clin-seq model: $clinseq_model_id");
    $clinseq_inputs_ok = 0;
  }
  
  return($clinseq_inputs_ok);
}


#Check whether instrument data of a particular type is available in a list of samples...
sub check_instrument_data{
  my $self = shift;
  my %args = @_;
  my $data_type = $args{'-data_type'};
  my @samples = @{$args{'-samples'}};
  
  my $samples_with_wgs = 0;
  my $samples_with_exome = 0;
  my $samples_with_normal_rnaseq = 0;
  my $samples_with_tumor_rnaseq = 0;

  my $normal_sample_common_names = $self->normal_sample_common_names;
  my $tumor_sample_common_names = $self->tumor_sample_common_names;

  foreach my $sample (@samples){
    my @normal_rnaseq;
    my @tumor_rnaseq;
    my @exome;
    my @wgs;
    my @unknown;
    my @other;
    my @sample_instrument_data = $sample->instrument_data;
    @sample_instrument_data = @{$self->exclude_instrument_data('-instrument_data'=>\@sample_instrument_data)};
    
    foreach my $instrument_data (@sample_instrument_data){
      my $trsn = $instrument_data->target_region_set_name;

      if ($sample->is_rna){
        my $scn = $sample->common_name || "NULL";
        push (@normal_rnaseq, $instrument_data) if ($scn =~ /$normal_sample_common_names/i);
        push (@tumor_rnaseq, $instrument_data) if ($scn =~ /$tumor_sample_common_names/i);
      }elsif ($trsn){
        my $fl = Genome::FeatureList->get(name => $trsn);
        if (not $fl or not $fl->content_type) {
          push @unknown, $instrument_data;
        }elsif ($fl->content_type eq 'exome') {
          push @exome, $instrument_data;
        }else {
          push @other, $instrument_data;
        }
      }else{
        push @wgs, $instrument_data;
      }
    }
    $samples_with_wgs++ if (scalar(@wgs));
    $samples_with_exome++ if (scalar(@exome));
    $samples_with_normal_rnaseq++ if (scalar(@normal_rnaseq));
    $samples_with_tumor_rnaseq++ if (scalar(@tumor_rnaseq));
  }

  my $samples_with_matching_data = 0;
  if ($data_type eq 'wgs'){
    $samples_with_matching_data = $samples_with_wgs;
  }elsif($data_type eq 'exome'){
    $samples_with_matching_data = $samples_with_exome;
  }elsif($data_type eq 'normal_rnaseq'){
    $samples_with_matching_data = $samples_with_normal_rnaseq;    
  }elsif($data_type eq 'tumor_rnaseq'){
    $samples_with_matching_data = $samples_with_tumor_rnaseq;        
  }else{
    $self->error_message("Data type specified to check_instrument_data not understood");
    exit 1;
  }

  return $samples_with_matching_data;
}


#Take an array of instrument data objects and return an amended array that removes certain instrument data if specified by the user with the --instrument_data_to_exclude parameter
sub exclude_instrument_data{
  my $self = shift;
  my %args = @_;
  my @instrument_data = @{$args{'-instrument_data'}};

  #If the user specified some instrument data to skip, deal with that first
  my @tmp1;
  if ($self->instrument_data_to_exclude){
    #Check format of list of data to exclude
    my @exclude_list = split(",", $self->instrument_data_to_exclude);
    unless (scalar(@exclude_list)){
      $self->error_message("Could not obtain instrument data IDs from list supplied by --instrument_data_to_exclude");
      exit 1;
    }
    #Cross reference instrument data to exclude with the instrument data supplied to the subroutine and create an amended array
    my %exclude_list;
    foreach my $iid (@exclude_list){
      $exclude_list{$iid}=1;
    }
    foreach my $instrument_data (@instrument_data){
      my $instrument_data_id = $instrument_data->id;
      push (@tmp1, $instrument_data) unless $exclude_list{$instrument_data_id};
    }
    @instrument_data = @tmp1;
  }

  #Now skip instrument data that is not really Illumina instrument data
  #TODO:  The following method may miss 'Solexa' instrument data that is imported but is not properly classed as 'Solexa' ...
  my @tmp2;
  foreach my $instrument_data (@instrument_data){
    next unless ($instrument_data->class eq "Genome::InstrumentData::Solexa");

    #TODO: Problem with the approach ... Some Solexa data produced here has been duplicated as 'imported data'
    #Its hard to tell the difference between this and real imported data that was not produced here. :(
    #Allow all instrument data that is defined as sequencing_platform of 'solexa' (includes those with class Genome::InstrumentData::Solexa and Genome::InstrumentData::Imported)
    #Limit to file types of bam or %fastq% (illumina fastq, sanger fastq, solexa fastq)
    #next unless ($instrument_data->sequencing_platform eq "solexa");
    #my $is_bam = 0;
    #my $is_fastq = 0;
    #if ($instrument_data->can("import_format")){
    #  $is_bam = 1 if ($instrument_data->import_format eq "bam");
    #  $is_fastq = 1 if ($instrument_data->import_format =~ /fastq/);
    #  next unless ($is_bam || $is_fastq);
    #}

    push(@tmp2, $instrument_data);
  }
  @instrument_data = @tmp2;

  #Now skip instrument data that has an index defined but where that index is defined as 'unknown'
  my @tmp3;
  foreach my $instrument_data (@instrument_data){
    my $index = $instrument_data->index_sequence;
    if ($index){
      next if ($index =~ /unknown/i);
    }
    push(@tmp3, $instrument_data);
  }
  @instrument_data = @tmp3;

  #Return the amended array
  return (\@instrument_data);
}


1;

