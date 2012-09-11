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

  #Are there suitable WGS reference alignment models in existence (If not, create)?  If so, what is the status?
  #- Are there tumor/normal DNA samples?
  #- Is there tumor/normal WGS instrument data? If so, is it all included in the existing model?


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


  $self->status_message("\n\n");

  return 1;
}


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
  $self->status_message("\nINDIVIDUAL");
  $self->status_message("id\tname\tcommon_name\tspecies_name\tupn\tgender\tethnicity\trace");
  $self->status_message("$id\t$name\t$common_name\t$species_name\t$upn\t$gender\t$ethnicity\t$race");
  return;
}


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
  $self->status_message("\nSAMPLES");
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


sub display_processing_profiles{
  my $self = shift;
  $self->status_message("\nDESIRED PROCESSING PROFILES");
  $self->status_message("reference-alignment: " . $self->ref_align_pp->name . " <" . $self->ref_align_pp->id . ">");
  $self->status_message("wgs somatic-variation: " . $self->wgs_somatic_variation_pp->name . " <" . $self->wgs_somatic_variation_pp->id . ">");
  $self->status_message("exome somatic-variation: " . $self->exome_somatic_variation_pp->name . " <" . $self->exome_somatic_variation_pp->id . ">");
  $self->status_message("rna-seq: " . $self->rnaseq_pp->name . " <" . $self->rnaseq_pp->id . ">");
  $self->status_message("clin-seq: " . $self->clinseq_pp->name . " <" . $self->clinseq_pp->id . ">");
  return;
}


sub display_inputs{
  my $self = shift;
  $self->status_message("\nDESIRED MODEL INPUTS");
  $self->status_message("reference_sequence_build: " . $self->reference_sequence_build->name);
  $self->status_message("reference_annotation_build: " . $self->annotation_build->name);
  $self->status_message("dbsnp_build: " . $self->dbsnp_build->__display_name__);
  $self->status_message("previously_discovered_variants: " . $self->previously_discovered_variations->__display_name__);
  return;
}



1;

