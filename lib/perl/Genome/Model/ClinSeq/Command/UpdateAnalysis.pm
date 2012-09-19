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
              require_user_verify => 0,
              doc => 'Sample(s) to target for clinseq analysis'
        },
        sample_type_filter => {
              is => 'Text',
              default => 'pcr product,pooled library',
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
              default => '127786607',
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
genome model clin-seq update-analysis --outdir='/tmp/update_analysis/' --individual='common_name=AML103'


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

  #Are there suitable WGS reference alignment models in existence (If not, create)?  If so, what is the status?
  #- Are there tumor/normal DNA samples?
  #- Is there tumor/normal WGS instrument data? If so, is it all included in the existing model?
  my @normal_wgs_ref_align_models;
  my @tumor_wgs_ref_align_models;
  my @wgs_somatic_variation_models;
  if (scalar(@dna_samples) == 2){
    $self->status_message("\nWGS REFERENCE-ALIGNMENT MODELS");
    @normal_wgs_ref_align_models = $self->check_ref_align_models('-data_type'=>'wgs', '-tissue_type'=>'normal', '-dna_samples'=>\@dna_samples);
    @tumor_wgs_ref_align_models = $self->check_ref_align_models('-data_type'=>'wgs', '-tissue_type'=>'tumor|met', '-dna_samples'=>\@dna_samples);
  
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
  my @normal_exome_ref_align_models;
  my @tumor_exome_ref_align_models;
  my @exome_somatic_variation_models;
  if (scalar(@dna_samples) == 2){
    $self->status_message("\nEXOME REFERENCE-ALIGNMENT MODELS");
    @normal_exome_ref_align_models = $self->check_ref_align_models('-data_type'=>'exome', '-tissue_type'=>'normal', '-dna_samples'=>\@dna_samples);
    @tumor_exome_ref_align_models = $self->check_ref_align_models('-data_type'=>'exome', '-tissue_type'=>'tumor|met', '-dna_samples'=>\@dna_samples);

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
    @normal_rnaseq_models = $self->check_rnaseq_models('-tissue_type'=>'normal', '-rna_samples'=>\@rna_samples);
    @tumor_rnaseq_models = $self->check_rnaseq_models('-tissue_type'=>'tumor', '-rna_samples'=>\@rna_samples);
  }

  #Is there a suitable clin-seq model in existence (If not, create)?  If so, what is the status?
  #- Only proceed with this is the prerequisite models (where data exists) are ready.


  #TODO: Remind the user whether this clinseq model is currently set as do-not-archive and remind them how to do that...

  $self->status_message("\n\n");

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
  
  my $sample_count = scalar(@samples);
  $self->status_message("\nEXAMINE SAMPLES");
  $self->status_message("Found $sample_count samples:");
  $self->status_message("id\tname\tcommon_name\tsample_type\tcell_type\ttissue_desc\tdefault_genotype_data_id\tmodel_count\tlibrary_count\tid_count");
  my $skip_count = 0;
  my $sample_mismatch = 0;
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
      $sample_mismatch++;
    }

    $self->status_message("$id\t$name\t$common_name\t$sample_type\t$cell_type\t$tissue_desc\t$default_genotype_data_id\t$model_count\t$library_count\t$id_count");
    push(@final_samples, $sample);
  }

  if ($sample_mismatch){
    $self->warning_message("\nFound $sample_mismatch samples provided by the user that do not match the specified patient.  Aborting ...\n");
    exit(1);
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


#Test for existence of DNA samples and return array of sample objects if found
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
  return (@dna_samples);
}


#Check for existence of RNA samples and return array of sample objects if found
sub rna_samples{
  my $self = shift;
  my %args = @_;
  my @samples = @{$args{'-samples'}};

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
  return (@rna_samples);
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
  my %trsns;

  foreach my $instrument_data (@sample_instrument_data){
    next unless ($instrument_data->class eq "Genome::InstrumentData::Solexa");
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
    $self->warning_message("Intrument data from more than one target region set are being combined...");
  }elsif($trsn_count == 0){
    $self->error_message("There is no instrument data with a target region set name!  How is this an exome data set?");
    exit(1);
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
    exit(1);
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
  }
  return $roi_name;
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
    $self->warning_message("\nModel: " . $model->id . " appears to be missing the following instrument data: @missing_data");
    $self->status_message("You should consider performing the following update before proceeding:");
    $self->status_message("genome model instrument-data assign --instrument-data='$id_string'");
    return 0;
  }
  return 1;
}


#Obtain a genotype microarray model object (if it exists) for a particular sample
sub get_genotype_microarray_model_id{
  my $self = shift;
  my %args = @_;
  my $sample = $args{'-sample'};

  #Get the imported data that should have been used for the genotype microarray model
  my $default_genotype_data = $sample->default_genotype_data;
  my $default_genotype_data_id = $default_genotype_data->id;
  
  my $genotype_microarray_model_id = 0;

  my @models = $sample->models;

  foreach my $model (@models){
    next unless ($model->class eq "Genome::Model::GenotypeMicroarray");
    $genotype_microarray_model_id = $model->id;
  }

  return $genotype_microarray_model_id;
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
    $self->error_message("\nDid not find a matching DNA sample for tissue type: $tissue_type");
  }elsif ($match > 1){
    $self->error_message("\nFound more than one matching DNA sample of tissue type: $tissue_type");
  }else{
    $self->status_message("\nFound a DNA sample " . $sample->name . " ($scn) matching tissue type: $tissue_type");
  }

  #Is there WGS or Exome data?  Get the instrument data of each type for this sample.
  #- return if there is no data of the desired type
  my @tmp;
  my @sample_instrument_data = $self->get_instrument_data('-sample'=>$sample, '-data_type'=>$data_type);
  return @tmp unless (scalar(@sample_instrument_data));

  my $subject_id = $sample->patient->id;
  my @models = $sample->models;
  my $model_count = scalar(@models);
  $self->status_message("\tStarting with " . $model_count . " $data_type models for this subject. Candidates that meet criteria:");

  #Test for correct processing profile, reference sequence build, annotation build, and dbsnp build
  #Also make sure that all the instrument data of wgs or exome type is being used (exome can be exome+wgs lanes)
  foreach my $model (@models){
    next unless ($model->class eq "Genome::Model::ReferenceAlignment");
    next unless ($model->processing_profile_id == $self->ref_align_pp->id);
    next unless ($model->reference_sequence_build->id == $self->reference_sequence_build->id);
    if ($model->can("annotation_reference_build")){
      next unless ($model->annotation_reference_build->id == $self->annotation_build->id);
    }else{
      next;
    }
    if ($model->can("dbsnp_build")){
      next unless ($model->dbsnp_build->id == $self->dbsnp_build->id);
    }else{
      next;
    }
  
    #Is this a WGS or an Exome model?
    #WGS models do not have exome data, Exome models have at least one lane of exome data - skip those that are not the current type being considered
    next unless ($self->determine_model_data_type('-model'=>$model) eq $data_type);

    #If the desired $data_type is exome.  Check that the TRSN and ROI have been set correctly, exclude models that are not
    next unless ($self->check_model_trsn_and_roi('-model'=>$model, '-data_type'=>$data_type));

    #Make sure all the wgs or exome data is associated with the model
    #In both wgs and exome models, additional data will be allowed to handle weird situations.  
    #Eventually we will need the ability to exclude data as well...
    next unless $self->check_for_missing_data('-model'=>$model, '-sample_instrument_data'=>\@sample_instrument_data);

    $self->status_message("\t\tName: " . $model->name);
    push (@final_models, $model);
  }

  my $final_model_count = scalar(@final_models);
  $self->status_message("\tFound " . $final_model_count . " suitable $data_type models (matching default or user specified criteria)");

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
    if ($tissue_type =~ /$current_scn/i){
      $match++;
      $sample = $s;
      $scn = $current_scn;
    }
  }
  if ($match == 0){
    $self->warning_message("\nDid not find a matching RNA sample for tissue type: $tissue_type");
    return @tmp;
  }elsif ($match > 1){
    $self->error_message("\nFound more than one matching RNA sample of tissue type: $tissue_type");
  }else{
    $self->status_message("\nFound an RNA sample " . $sample->name . " ($scn) matching tissue type: $tissue_type");
  }

  #Is there actually any RNA?
  #- return if there is no data of the desired type
  my @test = $sample->instrument_data;
  my @sample_instrument_data;
  foreach my $instrument_data (@test){
    next unless ($instrument_data->class eq "Genome::InstrumentData::Solexa");
    push (@sample_instrument_data, $instrument_data);
  }
  return @tmp unless (scalar(@sample_instrument_data));

  my $subject_id = $sample->patient->id;
  my @models = $sample->models;
  my $model_count = scalar(@models);
  $self->status_message("\tStarting with " . $model_count . " $tissue_type models for this subject. Candidates that meet criteria:");

  #Test for correct processing profile, reference sequence build, annotation build
  foreach my $model (@models){
    next unless ($model->class eq "Genome::Model::RnaSeq");
    next unless ($model->processing_profile_id == $self->rnaseq_pp->id);
    next unless ($model->reference_sequence_build->id == $self->reference_sequence_build->id);
    next unless ($model->annotation_build->id == $self->annotation_build->id);
  
    #Make sure all the rna-seq data is associated with the model
    next unless $self->check_for_missing_data('-model'=>$model, '-sample_instrument_data'=>\@sample_instrument_data);

    $self->status_message("\t\tName: " . $model->name);
    push (@final_models, $model);
  }

  my $final_model_count = scalar(@final_models);
  $self->status_message("\tFound " . $final_model_count . " suitable $tissue_type rna-seq models (matching default or user specified criteria)");

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
    exit(1);
  }

  #Get models for either tumor or normal DNA sample 
  my @existing_models;
  foreach my $sample(@dna_samples){
    my @models = $sample->models;
    push(@existing_models, @models);
  }
  my $existing_model_count = scalar(@existing_models);
  $self->status_message("\tStarting with " . $existing_model_count . " models for this pair of subjects. Candidates that meet criteria:");

  #Test for correct processing profile, annotation build, previously discovered variations
  foreach my $model (@existing_models){
    next unless ($model->class eq "Genome::Model::SomaticVariation");
    next unless ($model->processing_profile_id == $somatic_variation_pp_id);
    next unless ($model->annotation_build->id == $self->annotation_build->id);
    if ($model->can("previously_discovered_variations_build")){
      next unless ($model->previously_discovered_variations_build->id == $self->previously_discovered_variations->id);
    }else{
      next;
    }
    
    #Make sure one of the passing normal AND tumor reference alignment models are specified as inputs to the somatic variation model
    #next unless $self->check_for_missing_data('-model'=>$model, '-sample_instrument_data'=>\@sample_instrument_data);
    my $tumor_model_id = $model->tumor_model->id;
    my $tumor_model_match = 0;
    foreach my $model (@tumor_ref_align_models){
      $tumor_model_match = 1 if ($model->id == $tumor_model_id);
    }
    my $normal_model_id = $model->normal_model->id;
    my $normal_model_match = 0;
    foreach my $model (@normal_ref_align_models){
      $normal_model_match = 1 if ($model->id == $normal_model_id);
    }
    next unless ($tumor_model_match && $normal_model_match);

    $self->status_message("\t\tName: " . $model->name);
    push (@final_models, $model);
  }
  my $final_model_count = scalar(@final_models);
  $self->status_message("\tFound " . $final_model_count . " suitable $data_type somatic-variation models (matching default or user specified criteria)");

  #If there are no suitable models, one will need to be created
  unless ($final_model_count > 0){
    $self->create_somatic_variation_model('-processing_profile_id'=>$somatic_variation_pp_id, '-normal_models'=>\@normal_ref_align_models, '-tumor_models'=>\@tumor_ref_align_models);
    return @tmp;
  }

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

  #Get microarray model
  my $genotype_microarray_model_id = $self->get_genotype_microarray_model_id('-sample'=>$sample);

  #Come up with a descriptive model name.  May be more trouble than it is worth, perhaps it should be autogenerated?
  my $sample_name = $sample->name;
  my $annotation_id = $self->annotation_build->id;
  my $reference_build_id = $self->reference_sequence_build->id;
  my $ref_align_pp_id = $self->ref_align_pp->id;
  my $dbsnp_build_id = $self->dbsnp_build->id;

  my @commands;

  #WGS example
  if ($data_type eq 'wgs'){
    push(@commands, "\n#Create a WGS reference-alignment model as follows:");
    push(@commands, "genome model define reference-alignment  --reference-sequence-build='$reference_build_id'  --annotation-reference-build='$annotation_id'  --subject='$sample_name'  --processing-profile='$ref_align_pp_id'  --genotype-microarray-model='$genotype_microarray_model_id'  --dbsnp-build='$dbsnp_build_id'");
    push(@commands, "genome model instrument-data assign  --instrument-data='$iids_list'  --model=''");
    push(@commands, "genome model build start ''");
  }

  #Exome example
  if ($data_type eq 'exome'){
    #Determine the target region set name associated with this list of instrument-data
    my $target_region_set_name = $self->get_trsn('-instrument_data'=>\@sample_instrument_data);

    #Determine the desired region of interest set name for this target region set name
    my $region_of_interest_set_name = $self->get_roi_name('-target_region_set_name'=>$target_region_set_name);

    push(@commands, "\n#Create an Exome reference-alignment model as follows:");
    push(@commands, "genome model define reference-alignment  --reference-sequence-build='$reference_build_id'  --annotation-reference-build='$annotation_id'  --subject='$sample_name'  --processing-profile='$ref_align_pp_id'  --genotype-microarray-model='$genotype_microarray_model_id'  --dbsnp-build='$dbsnp_build_id'  --target-region-set-names='$target_region_set_name'  --region-of-interest-set-name='$region_of_interest_set_name'");
    push(@commands, "genome model instrument-data assign  --instrument-data='$iids_list'  --model=''");
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

  #Come up with a descriptive model name.  May be more trouble than it is worth, perhaps it should be autogenerated?
  my $sample_name = $sample->name;
  my $annotation_id = $self->annotation_build->id;
  my $reference_build_id = $self->reference_sequence_build->id;
  my $rnaseq_pp_id = $self->rnaseq_pp->id;

  my @commands;

  push(@commands, "\n#Create an RNA-seq model as follows:");
  push(@commands, "genome model define rna-seq  --reference-sequence-build='$reference_build_id'  --annotation-build='$annotation_id'  --subject='$sample_name'  --processing-profile='$rnaseq_pp_id'");
  push(@commands, "genome model instrument-data assign  --instrument-data='$iids_list'  --model=''");
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
  my $somatic_variation_pp_id = $self->somatic_variation_pp->id;
  my $annotation_build_id = $self->annotation_build->id;

  my @commands;
  push(@commands, "\n#Create a Somatic-Variation model as follows:");
  push(@commands, "genome model define somatic-variation  --processing-profile=$processing_profile_id  --tumor-model=$tumor_model_id  --normal-model=$normal_model_id  --annotation-build=$annotation_build_id");
  push(@commands, "genome model build start ''");

  return;
}


#Select the best model from a list of suitable models according to simple criteria that may depend on the model type being considered
sub select_best_model{
  my $self = shift;
  my %args = @_;
  my @models = @{'-models'};

  #If there is only one model, it is the best by default
  my $best_model = $models[0];
  return $best_model if (scalar(@models) == 1);

  #Rank models according to amount of instrument data available

  #If there is no clear winner, rank according to 'date_completed' of the last successful build?

  #If there is still no clear winner, arbitrarily return the first one in the array


  return $best_model;
}


#Check model build status for each model.  Return true if it is safe to proceed, launch builds where needed, report if we need to wait for running builds
sub check_models_status{
  my $self = shift;
  my %args = @_;
  my @models = @{$args{'-models'}};

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
        #Before attempting to use a model/build as input to another model, make sure it is not archived (and check whether it has been set as do-not-archive)
        #- Do not start automatically marking these things as do-not-archive, leave it to the end user to mark the actual desired clin-seq models as do-not-archive
        if ($build->is_archived){
          $self->status_message("Successful build $build_id of model $model_id that meets desired criteria is currently archived! Consider running: genome model build unarchive $build_id");
        }elsif ($build->archivable){
          $self->status_message("Successful build $build_id of model $model_id that meets desired criteria is currently archivable.  Consider running: genome model build set-do-not-archive $build_id");
        }
      }elsif ($status =~ /running/i){
        $running_builds++;
      }elsif ($status =~ /unstartable/i){
        $unstartable_builds++;
      }elsif ($status =~ /scheduled/i){
        $scheduled_builds++;
      }
    }

    #If a model has unstartable builds and now running or succeeded builds, warn the user
    if ($unstartable_builds && $running_builds == 0 && $succeeded_builds == 0){
      $self->status_message("\n\tWARNING\n\tModel: $model_id has unstartable builds only! Please investigate...\n\tgenome model status $model_id");
    }
    
    if ($running_builds){
      $self->status_message("\n\tWARNING\n\tModel: $model_id has running builds...\n\tgenome model status $model_id");
    }elsif ($succeeded_builds){
      $ready_model_count++;
    }elsif ($scheduled_builds){
      $self->status_message("\n\tWARNING\n\tModel: $model_id has a scheduled build...\n\tgenome model status $model_id");
    }else{
      $self->status_message("\n\tWARNING\n\tModel: $model_id needs a build ...\n\tgenome model status $model_id\n\tgenome model build start $model_id");
    }
  }

  #If all models are 'ready' return true, or perhaps if at least one model is 'ready' return true
  if ($ready_model_count > 0 && $ready_model_count == $model_count){
    return 1;
  }
  return undef;
}



1;

