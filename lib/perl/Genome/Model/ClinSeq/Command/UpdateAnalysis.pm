package Genome::Model::ClinSeq::Command::UpdateAnalysis;

#Written by Malachi Griffith

use strict;
use warnings;
use Genome;
use Data::Dumper;
use Term::ANSIColor qw(:constants);

class Genome::Model::ClinSeq::Command::UpdateAnalysis {
    is => 'Command::V2',
    has_input => [
        individual => { 
              is => 'Genome::Individual',
              is_many => 0,
              require_user_verify => 0,
              doc => 'Individual to query',
        },
        outdir => { 
              is => 'FilesystemPath',
              doc => 'Path where summary and commands will be written', 
        },
    ],
    has_optional => [
        samples => {
              is => 'Genome::Sample',
              is_many => 1,
              require_user_verify => 1,
              doc => 'Sample(s) to target for clinseq analysis'
        },
        sample_type_filter => {
              is => 'Text',
              default => 'pcr product',
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
              default => '2756469',
        },
        wgs_somatic_variation_pp => {
              is => 'Genome::ProcessingProfile',
              id_by => '_wgs_somatic_variation_pp_id',
              doc => 'Desired WGS somatic-variation processing profile',
        },
        _exome_somatic_variation_pp_id => {
              is => 'Number',
              default => '2756470',
        },
        exome_somatic_variation_pp => {
              is => 'Genome::ProcessingProfile',
              id_by => '_exome_somatic_variation_pp_id',
              doc => 'Desired Exome somatic-variation processing profile',
        },
        _rnaseq_pp_id => {
              is => 'Number',
              default => '2754795',
        },
        rnaseq_pp => {
              is => 'Genome::ProcessingProfile',
              id_by => '_rnaseq_pp_id',
              doc => 'Desired rna-seq processing profile',
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
              default => '106375969',
        },
        dbsnp_build => {
              is => 'Genome::Model::Build',
              id_by => '_dbsnp_build_id',
              doc => 'Desired dbSNP build',
        },
        _previously_discovered_variations_id => {
              is => 'Number',
              default => '110108854',
        },
        previously_discovered_variations => {
              is => 'Genome::Model::Build',
              id_by => '_previously_discovered_variations_id',
              doc => 'Desired previously discovered variants build',
        },
   ],
    doc => 'evaluate models/builds for an individual and help create/update a clinseq model that meets requested criteria',
};

sub help_synopsis {
    return <<EOS

When trying to determine what samples to use for the clin-seq analyis:

genome model clin-seq update-analysis --outdir='/tmp/update_analysis/' --individual='2878747495'
genome model clin-seq update-analysis --outdir='/tmp/update_analysis/' --individual='H_KA-306905'

genome model clin-seq update-analysis --outdir='/tmp/update_analysis/' --individual='H_KA-306905' --samples='id in [2878747496,2878747497,2879495575]'
genome model clin-seq update-analysis --outdir='/tmp/update_analysis/' --individual='H_KA-306905' --samples='name in ["H_KA-306905-1121472","H_KA-306905-1121474","H_KA-306905-S.4294"]'

All processing-profile and build input parameters can be specified by ID, name, etc. and should resolve

EOS
}

sub help_detail {
    return <<EOS
For a given individual, find the available samples and display to the user to select those desired for clinseq analysis

Once samples are set by the user, search for reference alignment, somatic variation, and rna-seq models for these samples

Determine whether these models meet certain criteria with respect to inputs and processing profile selections

EOS
}

sub __errors__ {
  my $self = shift;
  my @errors = $self->SUPER::__errors__(@_);

  unless (-e $self->outdir && -d $self->outdir) {
      push @errors, UR::Object::Tag->create(
	                                          type => 'error',
	                                          properties => ['outdir'],
	                                          desc => RED . "Outdir: " . $self->outdir . " not found or not a directory" . RESET,
                                          );
  }
  return @errors;
}

sub execute {
  my $self = shift;
  my $individual = $self->individual;
  my $outdir = $self->outdir;

  unless ($outdir =~ /\/$/){
    $outdir .= "/";
  }
  

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

  #Are there suitable WGS reference alignment models in existence (If not, create)?  If so, what is the status?
  #- Are there tumor/normal DNA samples?
  #- Is there tumor/normal WGS instrument data? If so, is it all included in the existing model?
  my @dna_samples = $self->dna_samples('-samples'=>\@samples);
  my @normal_wgs_ref_align_builds;
  my @tumor_wgs_ref_align_builds;
  if (scalar(@dna_samples)){
    @normal_wgs_ref_align_builds = $self->check_ref_align('-data_type'=>'wgs', '-tissue_type'=>'normal', '-dna_samples'=>\@dna_samples);
    @tumor_wgs_ref_align_builds = $self->check_ref_align('-data_type'=>'wgs', '-tissue_type'=>'tumor|met', '-dna_samples'=>\@dna_samples);
  }

  #Are there suitable Exome reference alignment models in existence (If not, create)? If so, what is the status?
  #- Are there tumor/normal DNA samples?
  #- Is there tumor/normal Exome instrument data?
  #- What is the target region set name (TRSN) and region of interest (ROI)? If there is data, is it all included in the existing model?


  #Is there a suitable rna-seq model in existence (If not, create)?  If so, what is the status?  
  #- Are there tumor/normal RNA samples?
  #- Is there tumor/normal RNA-seq instrument data? If so, is it all included in the existing model?


  #Is there a suitable WGS somatic variation model in existence (If not, create)?  If so, what is the status?
  #- Only proceed with this if the prerequisite WGS tumor/normal ref-align models exist


  #Is there a suitable Exome somatic variation model in existence
  #- Only proceed with this if the prerequisite Exome tumor/normal ref-align models exist


  #Is there a suitable clin-seq mode in existence (If not, create)?  If so, what is the status?  Does it contain neccessary instrument data?
  #- Only proeed with this is the prerequisite models (where data exists) are ready.



  #EXAMPLE MODEL QUERIES
  #my @models = Genome::Model->get(processing_profile_id => 2635769, subject_id => 2882207510);
  #my @models = Genome::Model::ReferenceAlignment->get(processing_profile_id => 2635769, subject_id => 2882207510, reference_sequence_build_id => 106942997);
  #my @models = Genome::Model::ReferenceAlignment->get(processing_profile_id => 2635769, subject_id => 2882207510, reference_sequence_build_id => 106942997, annotation_reference_build_id => 124434505);
  #my @models = Genome::Model::ReferenceAlignment->get(processing_profile_id => 2635769, subject_id => 2882207510, reference_sequence_build_id => 106942997, annotation_reference_build_id => 124434505, dbsnp_build_id => 106375969);



  $self->status_message("\n\n");

  return 1;
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
  
  my $sample_count = scalar(@samples);
  $self->status_message("\nEXAMINE SAMPLES");
  $self->status_message("Found $sample_count samples:");
  $self->status_message("id\tname\tcommon_name\tsample_type\tcell_type\ttissue_desc\tdefault_genotype_data_id\tmodel_count\tlibrary_count\tid_count");
  my $skip_count = 0;
  foreach my $sample (@samples){
    my $id = $sample->id;
    my $name = $sample->name;
    my $common_name = $sample->common_name || "NULL";
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
    unless ($patient_id == $individual_id){
      $self->error_message("ID of individual supplied by user ($individual_id) does not match that associated with a sample ($patient_id)");
    }

    $self->status_message("$id\t$name\t$common_name\t$sample_type\t$cell_type\t$tissue_desc\t$default_genotype_data_id\t$model_count\t$library_count\t$id_count");
    push(@final_samples, $sample);
  }
  if ($self->sample_type_filter){
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
  $self->status_message("clin-seq: " . $self->clinseq_pp->name . " <" . $self->clinseq_pp->id . ">");
  return;
}


#Display the model inputs to be targeted
sub display_inputs{
  my $self = shift;
  $self->status_message("\nEXAMINE DESIRED MODEL INPUTS");
  $self->status_message("reference_sequence_build: " . $self->reference_sequence_build->__display_name__);
  $self->status_message("annotation_build: " . $self->annotation_build->name . " (" . $self->annotation_build->id . ")");
  $self->status_message("dbsnp_build: " . $self->dbsnp_build->__display_name__);
  $self->status_message("previously_discovered_variations: " . $self->previously_discovered_variations->__display_name__);
  return;
}


#Test for existence of DNA samples and return sample objects if found
sub dna_samples{
  my $self = shift;
  my %args = @_;
  my @samples = @{$args{'-samples'}};
    
  my @dna_samples;
  foreach my $sample (@samples){
    if ($sample->sample_type =~ /dna/i){
      push(@dna_samples, $sample);
    }
  }
  my $dna_sample_count = scalar(@dna_samples);
  $self->status_message("Found " . $dna_sample_count . " DNA samples");
  if ($dna_sample_count){
    return (@dna_samples);
  }else{
    return 0;
  }
}


#Gather builds for a single subject that meet all the specified criteria
sub check_ref_align{
  my $self = shift;
  my %args = @_;
  my $data_type = $args{'-data_type'};
  my $tissue_type = $args{'-tissue_type'};
  my @dna_samples = @{$args{'-dna_samples'}};

  my @final_models;
  my @final_builds;

  #Make sure there is only one sample of this sample type (i.e. sample->common_name)
  my $match = 0;
  my $sample;
  my $scn;
  foreach my $s (@dna_samples){
    $scn = $s->common_name || "NULL";
    if ($tissue_type =~ /$scn/i){
      $match++;
      $sample = $s;
    }
  }
  if ($match == 0){
    $self->error_message("Did not find a matching sample for tissue type: $tissue_type");
  }elsif ($match > 1){
    $self->error_message("Found more than one matching sample of tissue type: $tissue_type");
  }else{
    $self->status_message("Found a sample " . $sample->name . " ($scn) matching tissue type: $tissue_type");
  }

  #Is there WGS or Exome data?  Get the instrument data of each type for this sample.
  #- return if there is no data of the desired type
  my @sample_instrument_data = $self->get_instrument_data('-sample'=>$sample, '-data_type'=>$data_type);
  return @final_builds unless (scalar(@sample_instrument_data));

  my $subject_id = $sample->patient->id;
  my @models = $sample->models;
  my $model_count = scalar(@models);
  $self->status_message("\tStarting with " . $model_count . " models for this subject");

  #Test for correct processing profile, reference sequence build, annotation build, and dbsnp build
  #Also make sure that all the instrument data of wgs or exome type is being used (exome can be exome+wgs lanes)
  foreach my $model (@models){
    next unless ($model->processing_profile_id == $self->ref_align_pp->id);
    next unless ($model->reference_sequence_build->id == $self->reference_sequence_build->id);
    next unless ($model->annotation_reference_build->id == $self->annotation_build->id);
    next unless ($model->dbsnp_build->id == $self->dbsnp_build->id);
  
    #Is this a WGS or an Exome model?
    #WGS models do not have exome data, Exome models have at least one lane of exome data
    next unless ($self->determine_model_data_type('-model'=>$model) eq $data_type);

    #TODO: Make sure all the wgs or exome data is associated with the model
    #In both wgs and exome models, additional data will be allowed to handle weird situations.  
    #Eventually we will need the ability to exclude data as well...
    #If some data is missing, allow the model to pass, but warn the user that they should update it
    $self->check_for_missing_data('-model'=>$model, '-sample_instrument_data'=>\@sample_instrument_data);


    $self->status_message("\t\tName: " . $model->name);
    push (@final_models, $model);
  }

  my $final_model_count = scalar(@final_models);
  $self->status_message("\tFound " . $final_model_count . " suitable models (matching default or user specified criteria)");

  return @final_builds;
}


#Get all instrument data for the sample of the specified type (wgs or exome)
sub get_instrument_data{
  my $self = shift;
  my %args = @_;
  my $sample = $args{'-sample'};
  my $data_type = $args{'-data_type'};

  my @sample_instrument_data = $sample->instrument_data;
  my $instrument_data_count = scalar(@sample_instrument_data);

  my @exome;
  my @wgs;
  my @unknown;
  my @other;
  foreach my $instrument_data (@sample_instrument_data){
    next unless ($instrument_data->class eq "Genome::InstrumentData::Solexa");
    my $trsn = $instrument_data->target_region_set_name;
    if ($trsn){
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

  if ($data_type eq 'wgs'){
    $self->status_message("Could not find any wgs data") unless (scalar(@wgs));
    return @wgs;
  }elsif($data_type eq 'exome'){
    $self->status_message("Could not find any exome data") unless (scalar(@exome));
    return @exome;
  }else{
    $self->error_message("Data type not understood in get_instrument_data");
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


#Compare instrument data on a model to that available for the sample (exome or wgs) and warn if data is missing
sub check_for_missing_data{
  my $self = shift;
  my %args = @_;
  my $model = $args{'-model'};
  my @sample_instrument_data = @{$args{'-sample_instrument_data'}};
  my @model_instrument_data = $model->instrument_data;
  my @missing_data;
  foreach my $sample_instrument_data (@sample_instrument_data){
    next unless ($sample_instrument_data->class eq "Genome::InstrumentData::Solexa");
    my $sid = $sample_instrument_data->id;
    my $match = 0;
    foreach my $model_instrument_data (@model_instrument_data){
      my $mid = $model_instrument_data->id;
      $match = 1 if ($mid == $sid);
    }
    push(@missing_data, $sid) unless $match;
  }
  if (scalar(@missing_data)){
    my $id_string = join(",", @missing_data);
    $self->status_message("\nWARNING!: Model: " . $model->id . " appears to be missing the following instrument data: @missing_data");
    $self->status_message("WARNING!: You should consider performing the following update before proceeding:");
    $self->status_message("genome model instrument-data assign --instrument-data='$id_string'");
  }
  return;
}


1;

