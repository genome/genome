package Genome::Model::MetagenomicCompositionShotgun;

use strict;
use warnings;

use Genome;
use File::Basename;
use File::Path 'make_path';
use Sys::Hostname;
use Data::Dumper;

our $UNALIGNED_TEMPDIR = '/tmp'; #'/gscmnt/sata844/info/hmp-mgs-test-temp';

class Genome::Model::MetagenomicCompositionShotgun {
    is => 'Genome::ModelDeprecated',
    has_param => [
        contamination_screen_pp_id => {
            is => 'Integer',
            doc => 'processing profile id to use for contamination screen',
            is_optional=> 1,
        },
        metagenomic_alignment_pp_id => {
            is => 'Integer',
            doc => 'processing profile id to use for metagenomic alignment',
        },
        unaligned_metagenomic_alignment_pp_id => {
            is => 'Integer',
            doc => 'processing profile id to use for realignment of unaligned reads from first metagenomic alignment',
            is_optional => 1,
        },
        first_viral_verification_alignment_pp_id => {
            is => 'Integer',
            doc => 'processing profile id to use for first viral verification alignment',
            is_optional => 1,
        },
        second_viral_verification_alignment_pp_id => {
            is => 'Integer',
            doc => 'processing profile id to use for first viral verification alignment',
            is_optional => 1,
        },
        merging_strategy => {
            is => 'Text',
            valid_values => [qw/ best_hit bwa /],
            doc => 'strategy used to merge results from metagenomic alignments. valid values : best_hit',
        },
        dust_unaligned_reads => {
            is => 'Boolean',
            default_value => 1,
            doc => 'flag determining if dusting is performed on unaligned reads from contamination screen step',
        },
        n_removal_threshold => {
            is => 'Integer',
            default_value => 0,
            doc => "Reads with this amount of n's will be removed from unaligned reads from contamination screen step before before optional dusting",
        },
        non_n_base_threshold => {
            is => 'Int',
            doc => 'reads with less than this amount of non-n bases will be removed in post-processing',
            default => 0,
        },
        mismatch_cutoff => {
            is => 'Integer',
            default_value=> 0,
            doc => 'mismatch cutoff (including softclip) for post metagenomic alignment processing',
        },
        include_taxonomy_report => {
            is => 'Boolean',
            default_value => 1,
            doc => 'When this flag is enabled, the model will attempt to grab taxonomic data for the metagenomic reports and produce a combined refcov-taxonomic final report.  Otherwise, only refcov will be run on the final metagenomic bam',
        },
        skip_qc_on_untrimmed_reads => {
            is => 'Boolean',
            default_value => 0,
            doc => "If this flag is specified, QC report will skip metric on the human-free, untrimmed data.",
        },
    ],
    has => [
        contamination_screen_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_mutable => 1,
            is_optional => 1,
            is_many => 1,
            via => 'inputs',
            to => 'value',
            where => [name => 'contamination_screen_reference', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence'],
        },
        unaligned_metagenomic_alignment_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_mutable => 1,
            via => 'inputs',
            to => 'value',
            where => [name => 'unaligned_metagenomic_alignment_reference ', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence'],
        },
        first_viral_verification_alignment_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_mutable => 1,
            is_optional => 1,
            via => 'inputs',
            to => 'value',
            where => [name => 'first_viral_verification_alignment_reference ', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence'],
        },
        second_viral_verification_alignment_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_mutable => 1,
            is_optional => 1,
            via => 'inputs',
            to => 'value',
            where => [name => 'second_viral_verification_alignment_reference ', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence'],
        },
        metagenomic_references => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_mutable => 1,
            is_many => 1,
            via => 'inputs',
            to => 'value',
            where => [name => 'metagenomic_alignment_reference', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence'],
        },
        contamination_screen_alignment_model => {
            is => 'Genome::Model::ReferenceAlignment',
            is_optional => 1,
            via => 'from_model_links',
            to => 'from_model',
            where => [role => 'contamination_screen_alignment_model'],
        },
        metagenomic_alignment_models => {
            is => 'Genome::Model::ReferenceAlignment',
            is_many => 1,
            via => 'from_model_links', 
            to => 'from_model',
            where => [role => 'metagenomic_alignment_model'],
        },
        unaligned_metagenomic_alignment_model => {
            is => 'Genome::Model::ReferenceAlignment',
            is_optional => 1,
            via => 'from_model_links',
            to => 'from_model',
            where => [role => 'unaligned_metagenomic_alignment_model'],
        },
        first_viral_verification_alignment_model => {
            is => 'Genome::Model::ReferenceAlignment',
            is_optional => 1, 
            via => 'from_model_links',
            to => 'from_model',
            where => [role => 'first_viral_verification_alignment_model'],
        },
        second_viral_verification_alignment_model => {
            is => 'Genome::Model::ReferenceAlignment',
            is_optional => 1, 
            via => 'from_model_links',
            to => 'from_model',
            where => [role => 'second_viral_verification_alignment_model'],
        },
        contamination_screen_pp => {
            is => 'Genome::ProcessingProfile::ReferenceAlignment',
            id_by => 'contamination_screen_pp_id',
            doc => 'processing profile to use for contamination screen',
            is_optional => 1,
        },
        metagenomic_alignment_pp => {
            is => 'Genome::ProcessingProfile::ReferenceAlignment',
            id_by => 'metagenomic_alignment_pp_id',
            doc => 'processing profile to use for metagenomic alignment',
        },
        unaligned_metagenomic_alignment_pp => {
            is => 'Genome::ProcessingProfile::ReferenceAlignment',
            id_by => 'unaligned_metagenomic_alignment_pp_id',
            doc => 'processing profile to use for realignment of unaligned reads from first metagenomic alignment',
            is_optional => 1,
        },
        first_viral_verification_alignment_pp => {
            is => 'Genome::ProcessingProfile::ReferenceAlignment',
            id_by => 'first_viral_verification_alignment_pp_id',
            doc => 'processing profile to use for first viral verification alignment',
            is_optional => 1,
        },
        second_viral_verification_alignment_pp => {
            is => 'Genome::ProcessingProfile::ReferenceAlignment',
            id_by => 'second_viral_verification_alignment_pp_id',
            doc => 'processing profile to use for first viral verification alignment',
            is_optional => 1,
        },
        sequencing_platform => {
            doc => 'The sequencing platform from whence the model data was generated',
            calculate_from => ['contamination_screen_pp'],
            calculate => q|$DB::single = 1; $contamination_screen_pp->sequencing_platform;|,
        },
    ],
};

sub build_subclass_name {
    return 'metagenomic-composition-shotgun';
}

sub delete {
    my $self = shift;
    for my $sub_model ($self->contamination_screen_alignment_model, $self->metagenomic_alignment_models) {
        next unless $sub_model;
        $sub_model->delete;
    }
    return $self->SUPER::delete(@_);
}

sub create{
    my $class = shift;

    my %params = @_;
    my $self = $class->SUPER::create(%params);
    return unless $self;

    if($self->contamination_screen_pp) {
        my $contamination_screen_model = $self->_create_underlying_contamination_screen_model();
        unless ($contamination_screen_model) {
            $self->error_message("Error creating contamination screening model!");
            $self->delete;
            return;
        }
    }

    my @metagenomic_models = $self->_create_underlying_metagenomic_models();
    my @metagenomic_references = $self->metagenomic_references;
    unless (@metagenomic_models == @metagenomic_references) {
        $self->error_message("Error creating metagenomic models!" . scalar(@metagenomic_models) . " " . scalar(@metagenomic_references));
        $self->delete;
        return;
    }

    if ($self->unaligned_metagenomic_alignment_pp) {
        my $unaligned_metagenomic_alignment_model = $self->_createunaligned_metagenomic_alignment_model();
        unless ($unaligned_metagenomic_alignment_model) {
            $self->error_message("Error creating unaligned metagenomic alignment model!");
            $self->delete;
            return;
        }
    }

    if ($self->first_viral_verification_alignment_pp) {
        my $first_viral_verification_alignment_model = $self->_createfirst_viral_verification_alignment_model();
        unless ($first_viral_verification_alignment_model) {
            $self->error_message("Error creating first viral verification alignment model!");
            $self->delete;
            return;
        }
    }

    if ($self->second_viral_verification_alignment_pp) {
        my $second_viral_verification_alignment_model = $self->_createsecond_viral_verification_alignment_model();
        unless ($second_viral_verification_alignment_model) {
            $self->error_message("Error creating second viral verification alignment model!");
            $self->delete;
            return;
        }
    }

    return $self;
}

sub _create_underlying_contamination_screen_model {
    my $self = shift;
    return $self->_create_model_for_type("contamination_screen");
}

sub _create_underlying_metagenomic_models {
    my $self = shift;

    my @new_objects;
    my $metagenomic_counter = 0;
    for my $metagenomic_reference ($self->metagenomic_references){ 
        $metagenomic_counter++;

        my %metagenomic_alignment_model_params = (
            processing_profile => $self->metagenomic_alignment_pp,
            subject_name => $self->subject_name, 
            name => $self->name.".metagenomic alignment model $metagenomic_counter",
            reference_sequence_build => $metagenomic_reference
        );
        my $metagenomic_alignment_model = Genome::Model::ReferenceAlignment->create( %metagenomic_alignment_model_params );

        unless ($metagenomic_alignment_model){
            $self->error_message("Couldn't create metagenomic reference model with params ".join(", " , map {$_ ."=>". $metagenomic_alignment_model_params{$_}} keys %metagenomic_alignment_model_params) );
            for (@new_objects){
                $_->delete;
            }
            return;
        }

        if ($metagenomic_alignment_model->reference_sequence_build($metagenomic_reference)){
            $self->status_message("updated reference sequence build on metagenomic alignment model ".$metagenomic_alignment_model->name);
        }else{
            $self->error_message("failed to update reference sequence build on metagenomic alignment model ".$metagenomic_alignment_model->name);
            for (@new_objects){
                $_->delete;
            }
            return;
        }

        push @new_objects, $metagenomic_alignment_model;
        $self->add_from_model(from_model=>$metagenomic_alignment_model, role=>'metagenomic_alignment_model');
        $self->status_message("Created metagenomic alignment model ".$metagenomic_alignment_model->__display_name__);
    }

    return @new_objects;
}

sub _createunaligned_metagenomic_alignment_model {
    my $self = shift;
    return $self->_create_model_for_type("unaligned_metagenomic_alignment");
}

sub _createfirst_viral_verification_alignment_model {
    my $self = shift;
    return $self->_create_model_for_type("first_viral_verification_alignment");
}

sub _createsecond_viral_verification_alignment_model {
    my $self = shift;
    return $self->_create_model_for_type("second_viral_verification_alignment");
}

sub _create_model_for_type {
    my $self = shift;
    my $type = shift;

    #CREATE UNDERLYING REFERENCE ALIGNMENT MODELS
    my $pp_accessor = $type."_pp";
    my $reference_accessor = $type."_reference";
    my %model_params = (
        processing_profile => $self->$pp_accessor,
        subject_name => $self->subject_name,
        name => $self->name.".$type model",
        reference_sequence_build=>$self->$reference_accessor,
    );
    my $model = Genome::Model::ReferenceAlignment->create( %model_params );

    unless ($model){
        $self->error_message("Couldn't create contamination screen model with params ".join(", ", map {$_ ."=>". $model_params{$_}} keys %model_params) );
        return;
    }

    if ($model->reference_sequence_build($self->$reference_accessor)){
        $self->status_message("updated reference sequence build on $type model ".$model->name);
    }else{
        $self->error_message("failed to update reference sequence build on $type model ".$model->name);
        return;
    }

    $self->add_from_model(from_model=> $model, role=>$type.'_alignment_model');
    $self->status_message("Created $type model ".$model->__display_name__);

    return $model;
}

sub _resolve_resource_requirements_for_build {
    my ($self, $build) = @_;
    my @instrument_data = $build->instrument_data;
    my $gtmp = 30 + 5 * (1 + scalar(@instrument_data));
    return "-R 'rusage[gtmp=$gtmp:mem=16000]' -M 16000000";
}


sub _execute_build {
    my ($self, $build) = @_;

    my $model = $build->model;
    $self->status_message('Build '.$model->__display_name__);

    # temp hack for debugging
    my $log_model_name = $model->name;
    $self->status_message("Starting build for model $log_model_name");

    my $screen_model;
    my $screen_build;
    my @screened_assignments;

    if ($self->contamination_screen_pp){ 
        $screen_model = $model->contamination_screen_alignment_model;
        unless ($screen_model) {
            die $self->error_message("couldn't grab contamination screen underlying model!");
        }

        # ENSURE WE HAVE INSTRUMENT DATA
        my @assignments = $model->instrument_data_inputs;
        if (@assignments == 0) {
            die $self->error_message("NO INSTRUMENT DATA ASSIGNED!");
        }

        # ASSIGN ANY NEW INSTRUMENT DATA TO THE CONTAMINATION SCREENING MODEL
        for my $assignment (@assignments) {
            my $instrument_data = $assignment->value;
            my $screen_assignment = $screen_model->instrument_data_input(
                value_id => $instrument_data->id 
            );
            if ($screen_assignment) {
                $self->status_message("Instrument data " . $instrument_data->__display_name__ . " is already assigned to the screening model");
            }
            else {
                my $cmd = 'genome model instrument-data assign --model-id '.$screen_model->id.' --instrument-data-id '.$instrument_data->id;
                if ($assignment->filter_desc){
                    $cmd.=' --filter '.$assignment->filter_desc;
                }
                my $exit_code = system($cmd);
                unless ($exit_code == 0){
                    die $self->error_message("Failed to add instrument data ".$instrument_data->id." to contamination screen model!");
                }
                my $confirm_screen_assignment = UR::Context->current->reload("Genome::Model::Input", model_id=>$screen_model->id, value_id=>$instrument_data->id);
                if ($confirm_screen_assignment){
                    $self->status_message("Added instrument data ".$instrument_data->id." to contamination screen model");
                }else{
                    die $self->error_message("Can't find instrument data assignment for ".$instrument_data->id." after attempting to assign to contamination screen model");
                }

            }
        }

        # BUILD HUMAN CONTAMINATION SCREEN MODEL
        $self->status_message("Building contamination screen model if necessary");
        $screen_build = $self->build_if_necessary_and_wait($screen_model);
        return if not $screen_build;

        my ($prev_from_build) = grep {$_->id eq $screen_build->id} $build->from_builds();
        $build->add_from_build(from_build=>$screen_build, role=>'contamination_screen_alignment_build') unless $prev_from_build;;

    }else{
        $self->status_message("Skipping contamination screen for instrument data");
    }

    # POST-PROCESS THE UNALIGNED READS FROM THE CONTAMINATION SCREEN MODEL

    my @imported_instrument_data_for_metagenomic_models;
    if ($self->contamination_screen_pp){

        #if we contamination screened, we'll use the alignment results, extract unaligned reads to a fastq and then post-process, create new imported instrument data
        $self->status_message("Processing and importing instrument data for any new unaligned reads");
        @screened_assignments = $screen_model->instrument_data_inputs;
        my @post_processed_unaligned_reads;
        for my $assignment (@screened_assignments) {
            my @alignment_results = $screen_build->alignment_results_for_instrument_data($assignment->value);
            if (@alignment_results > 1) {
                $self->error_message( "Multiple alignment_results found for instrument data assignment: " . $assignment->__display_name__);
                return;
            }
            if (@alignment_results == 0) {
                $self->error_message( "No alignment_results found for instrument data assignment: " . $assignment->__display_name__);
                return;
            }
            $self->status_message("Processing instrument data assignment ".$assignment->__display_name__." for unaligned reads import");

            my $alignment_result = $alignment_results[0];
            my @post_processed_unaligned_reads_for_alignment_result = $self->_process_unaligned_reads($alignment_result);
            push @post_processed_unaligned_reads, \@post_processed_unaligned_reads_for_alignment_result
        }

        unless (@post_processed_unaligned_reads == @screened_assignments) {
            Carp::confess("The count of post-processed unaligned reads does not match the count of screened instrument data assignments.");
        }
        @imported_instrument_data_for_metagenomic_models = @post_processed_unaligned_reads;

    }else{
        #if skipping contamination_screen, we need to extract the originally assigned imported instrument data and process and reimport it for the metagenomic screen.
        #sra data is stored in fastq/sangerqual format, so these just need to be extracted, dusted, n-removed
        my @sra_assignments = $build->instrument_data_inputs;
        for my $assignment(@sra_assignments){
            my @post_processed_reads = $self->_process_sra_instrument_data($assignment->value);
            #TODO, this doesn't need to be an array of array refs, since we should get a 1 to 1 sra inst_data to post-processed imported inst_data, but this is how it's done in the original pipeline, where the 1 to 1 convention doesn't hold.  We're sticking to this standard for now.
            push @imported_instrument_data_for_metagenomic_models, \@post_processed_reads;
        }
        @screened_assignments = @sra_assignments;
    }


    # ASSIGN THE POST-PROCESSED READS TO THE METAGENOMIC MODELS
    my @metagenomic_models = $model->metagenomic_alignment_models;
    for my $metagenomic_model (@metagenomic_models) {
        my %assignments_expected;
        for my $n (0..$#imported_instrument_data_for_metagenomic_models) {
            my $prev_assignment = $screened_assignments[$n];
            my $filter_desc = $prev_assignment->filter_desc;
            my $post_processed_instdata_for_prev_assignment = $imported_instrument_data_for_metagenomic_models[$n];
            for my $instrument_data (@$post_processed_instdata_for_prev_assignment) {
                my $metagenomic_assignment = $metagenomic_model->instrument_data_input(
                    value_id => $instrument_data->id #TODO inst data switch
                );
                if ($metagenomic_assignment) {
                    $self->status_message("Instrument data " . $instrument_data->__display_name__ . " is already assigned to model " . $metagenomic_model->__display_name__);
                }
                else {
                    my $cmd = 'genome model instrument-data assign --model-id '.$metagenomic_model->id.' --instrument-data-id '.$instrument_data->id;
                    if ($filter_desc){
                        $cmd.=' --filter '.$filter_desc;
                    }
                    my $exit_code = system($cmd);
                    unless ($exit_code == 0){
                        die $self->error_message("Failed to add instrument data ".$instrument_data->id." to model ".$metagenomic_model->__display_name__."!");
                    }
                    $metagenomic_assignment = UR::Context->current->reload("Genome::Model::Input", model_id=>$metagenomic_model->id, value_id=>$instrument_data->id);
                    if ($metagenomic_assignment){
                        $self->status_message("Added instrument data ".$instrument_data->id." to model ".$metagenomic_model->__display_name__);
                    }else{
                        die $self->error_message("Can't find instrument data assignment for ".$instrument_data->id." after attempting to assign to model ".$metagenomic_model->__display_name__);
                    }
                }
                $assignments_expected{$metagenomic_assignment->id} = $metagenomic_assignment;
            }
        }

        #TODO  When a model with skip_contamination_screen restarts, the potential extra fragment instrument data from n-removal will not be returned from _process_sra_data, and therefore won't be an expected assignment, so we'll need to check if the extra assignment's source_data_files contains(is derived from) the original instrument data assigned to the model. tldr: don't break on skip_contamination_restarts if extra inst_data looks legit.
        my @mcs_instrument_data_ids = map {$_->id} $build->instrument_data;

        # ensure there are no other odd assignments on the model besides those expected
        # this can happen if instrument-data is re-processed (deleted)
        for my $assignment ($metagenomic_model->instrument_data_inputs) {
            unless ($assignments_expected{$assignment->id}) {
                my $instrument_data = $assignment->instrument_data; #TODO inst data switch to value
                if ($instrument_data) {
                    my ($derived_from) = $instrument_data->original_data_path =~ m{^/tmp/(?:unaligned_reads/)?(\d+)}; 
                    #if unaligned reads is in the data path, the id may be an instrument data id(deprecated way of storing original data path), or an alignment result id, here we will figure out which and return the derived from instrument data id
                    my $inst_data = Genome::InstrumentData->get($derived_from);
                    unless($inst_data){
                        my $alignment =Genome::InstrumentData::AlignmentResult->get($derived_from);
                        unless ($alignment){
                            die $self->error_message("Couldn't determine derived_from instrument data id from original data path ".$instrument_data->original_data_path);
                        }
                        $derived_from = $alignment->instrument_data->id;
                    }
                    unless (grep {$derived_from eq $_} @mcs_instrument_data_ids){
                        $self->error_message(
                            "Odd assignment found on model "
                            . $metagenomic_model->__display_name__
                            . " for instrument data "
                            . $instrument_data->__display_name__
                        );
                        Carp::confess($self->error_message);
                    }
                }
                else {
                    $self->warning_message(
                        "Odd assignment found on model "
                        . $model->__display_name__
                        . " for MISSING instrument data.  Deleting the assignment."
                    );
                    $assignment->delete;
                }
            }
        }
    }


    # BUILD THE METAGENOMIC REFERENCE ALIGNMENT MODELS
    my @metagenomic_builds = $self->build_if_necessary_and_wait(@metagenomic_models);
    return if not @metagenomic_builds;

    for my $meta_build (@metagenomic_builds){
        my ($prev_meta_from_build) = grep {$_->id eq $meta_build->id} $build->from_builds();
        $build->add_from_build(from_build=>$meta_build, role=>'metagenomic_alignment_build') unless $prev_meta_from_build;
    }
    # SYMLINK ALIGNMENT FILES TO BUILD DIRECTORY
    my $data_directory = $build->data_directory;
    if ($self->contamination_screen_pp_id){
        my ($screen_bam, $screen_flagstat) = $self->get_bam_and_flagstat_from_build($screen_build);

        unless ($screen_bam and $screen_flagstat){
            die $self->error_message("Bam or flagstat undefined for contamination screen build(screen bam: $screen_bam, screen flagstat: $screen_flagstat) ");
        }
        unless (-e $screen_bam and -e $screen_flagstat){
            die $self->error_message("Bam or flagstat doesn't exist for contamination screen build(screen bam: $screen_bam, screen flagstat: $screen_flagstat) ");
        }
        $self->symlink($screen_bam, "$data_directory/contamination_screen.bam");
        $self->symlink($screen_flagstat, "$data_directory/contamination_screen.bam.flagstat");
    }

    my $counter;
    my @meta_bams;
    for my $meta_build (@metagenomic_builds){
        $counter++;
        my ($meta_bam, $meta_flagstat) = $self->get_bam_and_flagstat_from_build($meta_build);
        push @meta_bams, $meta_bam;

        unless ($meta_bam and $meta_flagstat){ 
            die $self->error_message("Bam or flagstat undefined for metagenomic alignemnt build $counter(meta bam: $meta_bam, meta flagstat: $meta_flagstat) ");
        }
        unless (-e $meta_bam and -e $meta_flagstat){
            die $self->error_message("Bam or flagstat doesn't exist for metagenomic alignment build $counter(meta bam: $meta_bam, meta flagstat: $meta_flagstat) ");
        }

        $self->symlink($meta_bam, "$data_directory/metagenomic_alignment$counter.bam");
        $self->symlink($meta_flagstat, "$data_directory/metagenomic_alignment$counter.bam.flagstat");
    }

    # REPORTS

    # enable "verbose" logging so we can actually see status messages from these methods
    local $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;

    if ($self->contamination_screen_pp_id){
        my $qc_report = Genome::Model::MetagenomicCompositionShotgun::Command::QcReport->create(
            build_id => $build->id,
            skip_qc_on_untrimmed_reads => $self->skip_qc_on_untrimmed_reads,
        );
        unless($qc_report->execute()) {
            die $self->error_message("Failed to create QC report!");
        }
    }

    my $final_metagenomic_bam = $build->_final_metagenomic_bam();
    if (@meta_bams > 1){
        $final_metagenomic_bam = $self->merge_metagenomic_bams(\@meta_bams, $final_metagenomic_bam);
    }else{
        $self->symlink($meta_bams[0], $final_metagenomic_bam);
    }

    # TODO we should probably lock here for all instrument data id's involved
    my $locks_ref;

    # If we are trying to re-align unaligned stuff from the previous metagenomic alignment...
    # Get the unaligned reads from the metagenomic alignment, make a model, and re-align
    my @realigned_instrument_data;
    my %metagenomic_instrument_data_to_alignment;
    if ($self->unaligned_metagenomic_alignment_pp) {
        for my $metagenomic_alignment_build ($build->_metagenomic_alignment_builds){
            for my $alignment_result ($metagenomic_alignment_build->get_alignments){
                my $instrument_data_id = $alignment_result->instrument_data->id;
                $metagenomic_instrument_data_to_alignment{$instrument_data_id} = $alignment_result->id;
            }
        }

        #EXTRACT UNALIGNED READS FROM METAGENOMIC ALIGNMENT AND IMPORT NEW INSTRUMENT DATA(1)
        my ($orig_data_paths_to_unaligned_fastq_files) = $self->extract_unaligned_reads_from_bam($final_metagenomic_bam, \%metagenomic_instrument_data_to_alignment);
        # $original_data_path_to_fastq_files{$original_data_path}->{file} = $file;_
        # $original_data_path_to_fastq_files{$original_data_path}->{instrument_data_id} = $instrument_data_id;
        my @unaligned_instrument_data = $self->upload_instrument_data_and_unlock($orig_data_paths_to_unaligned_fastq_files, $locks_ref);

        #ASSIGN INSTRUMENT DATA TO UNALIGNED_METAGENOMIC_ALIGNMENT_MODEL
        my $unaligned_metagenomic_alignment_model = $model->unaligned_metagenomic_alignment_model;
        unless ($unaligned_metagenomic_alignment_model) {
            die $self->error_message("couldn't grab unaligned metagenomic alignment underlying model!");
        }
        ## Assign unaligned instrument data from (1) to this model 
        for my $instrument_data (@unaligned_instrument_data) {
            $unaligned_metagenomic_alignment_model->add_instrument_data(
                value => $instrument_data,
                #filter_desc => "What do we put here?",  FIXME
            );
        }

        my $unaligned_metagenomic_alignment_build = $self->build_if_necessary_and_wait($unaligned_metagenomic_alignment_model);
        return if not $unaligned_metagenomic_alignment_build;

        my $realigned_bam = $unaligned_metagenomic_alignment_build->whole_rmdup_bam_file;

        Genome::Sys->create_symlink($realigned_bam, $build->data_directory . "/realigned.bam");
        Genome::Sys->create_symlink($realigned_bam.".flagstat", $build->data_directory . "/realigned.bam.flagstat");

        #EXTRACT ALIGNED (AND VIRAL) READS FROM REALIGNING UNALIGNED_METAGENOMIC_ALIGNMENT_MODEL AND IMPORT(2)
        my %realigned_instrument_data_to_alignment;
        for my $alignment_result ($unaligned_metagenomic_alignment_build->alignment_results){
            my $instrument_data_id = $alignment_result->instrument_data->id;
            $realigned_instrument_data_to_alignment{$instrument_data_id} ||= [];
            push @{$realigned_instrument_data_to_alignment{$instrument_data_id}}, $alignment_result->id;
        }

        my ($orig_data_paths_to_realigned_fastq_files) = $self->extract_aligned_reads_from_bam($realigned_bam, \%realigned_instrument_data_to_alignment); # @fastq_files includes fragment, forward and reverse... which do we use?
        $orig_data_paths_to_realigned_fastq_files = $self->extract_virus_reads_from_fastqs($orig_data_paths_to_realigned_fastq_files);
        @realigned_instrument_data = $self->upload_instrument_data_and_unlock($orig_data_paths_to_realigned_fastq_files, $locks_ref);
    }

    # If we are doing viral realignment for either model...
    # Assign the previously aligned metagenomic reads (after screening for only viral reads) ... and the re-aligned reads if we did them
    if ($self->first_viral_verification_alignment_pp || $self->second_viral_verification_alignment_pp) {
        #EXTRACT ALIGNED (AND VIRAL) READS FROM THE ORIGINAL METAGENOMIC ALIGNMENT AND IMPORT(3)
        my ($orig_data_paths_to_metagenomic_fastq_files) = $self->extract_aligned_reads_from_bam($final_metagenomic_bam, \%metagenomic_instrument_data_to_alignment); 
        $orig_data_paths_to_metagenomic_fastq_files = $self->extract_virus_reads_from_fastqs($orig_data_paths_to_metagenomic_fastq_files);
        my @aligned_instrument_data = $self->upload_instrument_data_and_unlock($orig_data_paths_to_metagenomic_fastq_files, $locks_ref); 

        #ASSIGN NEW IMPORTED READS(2,3) TO VIRAL_VERIFICATION_ALIGNMENT MODELS AND BUILD
        if ($self->first_viral_verification_alignment_pp) {
            my $virus_alignment_model_1 = $model->first_viral_verification_alignment_model;
            unless ($virus_alignment_model_1) {
                die $self->error_message("couldn't grab virus alignment underlying model!");
            }

            for my $instrument_data (@realigned_instrument_data, @aligned_instrument_data) {
                $virus_alignment_model_1->add_instrument_data(
                    value => $instrument_data,
                    # filter_desc => "What do we put here?", FIXME
                );
            }
            my $virus_alignment_build_1 = $self->build_if_necessary_and_wait($virus_alignment_model_1);
            return if not $virus_alignment_build_1;
            my $virus_1_bam = $virus_alignment_build_1->whole_rmdup_bam_file;
            Genome::Sys->create_symlink($virus_1_bam, $build->data_directory . "/virus_1.bam");
            Genome::Sys->create_symlink($virus_1_bam.".flagstat", $build->data_directory . "/virus_1.bam.flagstat");
        }

        if ($self->second_viral_verification_alignment_pp) {
            my $virus_alignment_model_2 = $model->second_viral_verification_alignment_pp;
            unless ($virus_alignment_model_2) {
                die $self->error_message("couldn't grab virus alignment underlying model!");
            }
            for my $instrument_data (@realigned_instrument_data, @aligned_instrument_data) {
                $virus_alignment_model_2->add_instrument_data(
                    value => $instrument_data,
                    # filter_desc => "What do we put here?", FIXME
                );
            }
            my $virus_alignment_build_2 = $self->build_if_necessary_and_wait($virus_alignment_model_2);
            return if not $virus_alignment_build_2;
            my $virus_2_bam = $virus_alignment_build_2->whole_rmdup_bam_file;
            Genome::Sys->create_symlink($virus_2_bam, $build->data_directory . "/virus_2.bam");
            Genome::Sys->create_symlink($virus_2_bam.".flagstat", $build->data_directory . "/virus_2.bam.flagstat");
        }
    }

    # CREATE BAM INDEX
    my $bam_file = $build->_final_metagenomic_bam;
    my $samtools_version = $self->metagenomic_alignment_pp->samtools_version;
    my $cmd = Genome::Model::Tools::Sam::IndexBam->create(
        bam_file => $bam_file,
        use_version => $samtools_version,
    );

    unless($cmd->execute){
        $self->error_message("Failed to index bam file $bam_file !");
        return;
    }

    # REFCOV
    my $refcov = $self->_run_refcov($build);
    if ( not $refcov ) {
        $self->error_message('Failed to run refcov');
        return;
    }

    # TAXONOMY REPORT
    if ( $self->include_taxonomy_report ){
        my $taxonomy = $self->_run_taxonomy_report($build);
        if ( not $taxonomy ) {
            $self->error_message('Failed to run taxonomy report');
            return;
        }
    }

    # VALIDATE
    unless($self->contamination_screen_pp){ #TODO: update validate to deal with this arg
        my $validate_build = Genome::Model::MetagenomicCompositionShotgun::Command::Validate->create(build_id => $build->id);
        $validate_build->dump_status_messages(1);
        unless($validate_build->execute()) {
            die $self->error_message("Failed to validate build!");
        }
    }

    $self->status_message('Build success!');

    return 1;
}

sub _run_refcov {
    my ($self, $build) = @_;

    $self->status_message('Run Refcov');

    my $refcov_output = $build->refcov_output;
    $self->status_message('Refcov output: '.$refcov_output);
    if (-e $refcov_output and -e "$refcov_output.ok"){
        $self->status_message("Refcov already complete, shortcutting");
        return 1;
    }

    my $sorted_bam = $build->_final_metagenomic_bam;
    $self->status_message('Metagenomic BAM: '.$sorted_bam);
    if ( not -s $sorted_bam ) {
        die $self->error_message('Sorted BAM does not exist!');
    }

    my $regions_file = $build->metagenomic_reference_regions_file; 
    $self->status_message('Regions file: '.$regions_file);
    if ( not -s $regions_file ) {
        die $self->error_message("No regions file ($regions_file) does not exist");
    }

    my $command = "genome-perl5.10 -S gmt ref-cov standard";
    $command .= " --alignment-file-path ".$sorted_bam;
    $command .= " --roi-file-path ".$regions_file;
    $command .= " --stats-file ".$refcov_output;
    $command .= " --min-depth-filter 0";

    $self->status_message($command);
    my $rv = eval{
        Genome::Sys->shellcmd(
            cmd=>$command,
            #output_files=>[$refcov_output],
            #input_files => [$sorted_bam, $regions_file],
        );
    };
    if ( not $rv ) {
        $self->error_message("Failed to run refcov command: $@");
        return;
    }
    if ( not -e $refcov_output ) {
        $self->error_message("Refcov command succeeded, but output file ($refcov_output) does not exist");
        return;
    }

    $self->status_message('Run Refcov...OK');

    return 1;
}

sub _run_taxonomy_report {
    my ($self, $build) = @_;

    $self->status_message('Run taxonomy report');

    my $report = Genome::Model::MetagenomicCompositionShotgun::Command::TaxonomyReport->create(build => $build);
    if ( not $report ) {
        $self->error_message('Failed to create taxnomy report command');
        return;
    }
    $report->dump_status_messages(1);
    unless( $report->execute ) {
        $self->error_message('Failed to exectue taxonomy report command');
        return;
    }

    $self->status_message('Run taxonomy report...OK');

    return 1;
}

sub merge_metagenomic_bams{
    my ($self, $meta_bams, $sorted_bam) = @_;
    if (-e $sorted_bam and -e $sorted_bam.".OK"){
        $self->status_message("sorted metagenomic merged bam already produced, skipping");
    }else{
        my $merged_bam = $sorted_bam.".name_sorted.bam";
        $self->status_message("starting sort and merge");

        my $sort_and_merge_meta = Genome::Model::Tools::Sam::SortAndMergeSplitReferenceAlignments->create(
            input_files => $meta_bams,
            output_file => $merged_bam,
        );
        unless($sort_and_merge_meta->execute()) {
            die $self->error_message("Failed to sort and merge metagenomic bams: $@");
        }

        unless (-s $merged_bam){
            die $self->error_message("Merged bam has no size!");
        }

        $self->status_message("starting position sort of merged bam");

        my $sort_merged_bam = Genome::Model::Tools::Sam::SortBam->create(
            file_name => $merged_bam,
            output_file => $sorted_bam,
        );
        unless($sort_merged_bam->execute()) {
            die $self->error_message("Failed to position sort merged metagenomic bam.");
        }

        unless (-s $sorted_bam){
            die $self->error_message("Sorted bam has no size!");
        }

        system ("touch $sorted_bam.OK");
        unlink $merged_bam;
    }
    return $sorted_bam;
}

# TODO This method should take a list of fastq files and return new fastq files that contain only viral reads.
# The logic for this method is not yet available. So this serves as a placeholder.
sub extract_virus_reads_from_fastqs {
    my $self = shift;
    return @_;
}

sub extract_unaligned_reads_from_bam {
    my $self = shift;
    return $self->extract_reads_from_bam(@_, 0);
}

sub extract_aligned_reads_from_bam {
    my $self = shift;
    return $self->extract_reads_from_bam(@_, 1);
}

# Extracts reads from a bam file and checks the output
# $extract_aligned_reads should be provided a 0 to extract unaligned reads, or a 1 to extract aligned reads
sub extract_reads_from_bam {
    my ($self, $bam, $instrument_data_to_alignment, $extract_aligned_reads) = @_;
    my ($frag_fastq, $fwd_fastq, $rev_fastq);

    #TODO we should check to see if we have already done this work and uploaded instrument data before doing it again

    my $output_dir = Genome::Sys->create_temp_directory;

    unless (-s $bam) {
        die $self->error_message("Failed to find expected BAM file $bam\n");
    }

    $self->status_message("Preparing imported instrument data for import path $output_dir");
#    my $extract_unaligned = Genome::Model::Tools::BioSamtools::BamToUnalignedFastq->create(
#        bam_file => $bam,
#        output_directory => $output_dir,
#        print_aligned => $extract_aligned_reads,
#        );
#    $self->execute_or_die($extract_unaligned);

    my $cmd = "genome-perl5.10 -S gmt bio-samtools bam-to-unaligned-fastq --bam-file $bam --output-directory $output_dir";
    if($extract_aligned_reads) {
        $cmd .= " --print-aligned";
    }
    if(system($cmd)) {
        die("$cmd failed: $?");
    }

    my @files = glob("$output_dir/*/*_sequence.txt");


    my %original_data_path_to_fastq_files;

    my $dir = Genome::Sys->create_temp_directory();

    # Create a hash to map original_data_path -- (to be used when uploading , describing where the reads came from (which alignment result) and what they are, (aligned on unaligned) )
    # -- to the original location of the file, and the original instrument data id it came from
    for my $file (@files) {
        my ($instrument_data_id,$lane,$pe_segment) = $file =~ m{(\d+)/s_(\d)_?([12])?_sequence.txt};
        my $alignment_result = $instrument_data_to_alignment->{$instrument_data_id};
        my $read_type;
        if ($extract_aligned_reads) {
            $read_type = "aligned_reads";
        } else {
            $read_type = "unaligned_reads";
        }
        my $part = "$dir/$read_type/$alignment_result";
        make_path($part);
        if ($pe_segment){
            my $original_data_path = "$part/s_".$lane."_1_sequence.txt," . "$part/s_".$lane."_2_sequence.txt";
            $original_data_path_to_fastq_files{$original_data_path}->{file}->{$pe_segment} = $file;
            $original_data_path_to_fastq_files{$original_data_path}->{instrument_data_id} = $instrument_data_id;
        }else{
            my $original_data_path = "$part/s_".$lane."_sequence.txt";
            $original_data_path_to_fastq_files{$original_data_path}->{file} = $file;
            $original_data_path_to_fastq_files{$original_data_path}->{instrument_data_id} = $instrument_data_id;
        }
    }

    $self->status_message("Extracted unaligned reads from bam file($bam)");
    return \%original_data_path_to_fastq_files;
}

sub get_imported_instrument_data_or_upload_paths {  #TODO not finished, not currently used
    my ($self, $orig_inst_data, @paths) = @_;
    my @inst_data;
    my @upload_paths;
    my @locks;
    my $tmp_dir;
    my $subdir;

    # check for previous unaligned reads
    $self->status_message("Checking for previously imported unaligned and post-processed reads from: $tmp_dir/$subdir");
    for my $path (@paths) {
        my $inst_data = Genome::InstrumentData::Imported->get(original_data_path => $path);
        if ($inst_data) {
            $self->status_message("imported instrument data already found for path $path, skipping");
            push @inst_data, $inst_data;
        }
        else {
            for my $sub_path (split(',', $path)) {
                my $lock = $self->lock($orig_inst_data->id, basename($sub_path));
                if ($lock) {
                    $self->status_message("Locked $sub_path successfully.");
                } else {
                    die $self->error_message("Failed to lock $sub_path.");
                }
                push @locks, $lock;
            }
            push @upload_paths, $path;
        }
    }

    if (@upload_paths) {
        for my $path (@upload_paths) {
            $self->status_message("planning to upload for $path");
        }
        return @upload_paths;
    }
    else {
        $self->status_message("skipping read processing since all data is already processed and uploaded");
        return @inst_data;
    }
    return \@inst_data, \@upload_paths;
}

sub upload_instrument_data_and_unlock {
    my ($self, $orig_data_paths_to_fastq_files, $locks_ref) = @_;
    #TODO: Actually use locks.
    my @locks;# = @$locks_ref;

    my @properties_from_prior = qw/ run_name subset_name sequencing_platform median_insert_size sd_above_insert_size library_name sample_name /;
    my @instrument_data;
    for my $original_data_path (keys %$orig_data_paths_to_fastq_files) {
        # If it is paired end
        if ($original_data_path =~ /,/) {
            my ($original_data_path_1, $original_data_path_2) = split ",", $original_data_path;
            Genome::Sys->create_symlink($orig_data_paths_to_fastq_files->{$original_data_path}->{file}->{1}, $original_data_path_1);
            Genome::Sys->create_symlink($orig_data_paths_to_fastq_files->{$original_data_path}->{file}->{2}, $original_data_path_2);
            if(! -s $original_data_path_1) {
                $self->warning_message('$original_data_path_1 is empty');
                next;
            }
            if(! -s $original_data_path_2) {
                $self->warning_message('$original_data_path_2 is empty');
                next;
            }
        } else {
            Genome::Sys->create_symlink($orig_data_paths_to_fastq_files->{$original_data_path}->{file}, $original_data_path);
            if(! -s $original_data_path) {
                $self->warning_message('$original_data_path is empty');
                next;
            }
        }
        # Link the actual file to the descriptive location in original_data_path
        my $orig_inst_data = Genome::InstrumentData->get($orig_data_paths_to_fastq_files->{$original_data_path}->{instrument_data_id});

        $self->status_message("uploading new instrument data ...");
        my @errors;
        my %properties_from_prior;
        for my $property_name (@properties_from_prior) {
            my $value = $orig_inst_data->$property_name;
            no warnings;
            $self->status_message("Value for $property_name is $value");
            $properties_from_prior{$property_name} = $value;
        }

        if ($original_data_path =~ /,/){
            $properties_from_prior{is_paired_end} = 1;
        }else{
            $properties_from_prior{is_paired_end} = 0;
        }
        my %params = (
            %properties_from_prior,
            source_data_files => $original_data_path,
            import_format => 'illumina fastq',
        );
        $self->status_message("importing fastq with the following params:" . Data::Dumper::Dumper(\%params));


        my $command = Genome::InstrumentData::Command::Import::Fastq->create(%params);
        unless ($command) {
            $self->error_message( "Couldn't create command to import unaligned fastq instrument data!");
        };
        my $result = $command->execute();
        unless ($result) {
            die $self->error_message( "Error importing data from $original_data_path! " . Genome::InstrumentData::Command::Import::Fastq->error_message() );
        }
        $self->status_message("committing newly created imported instrument data");
        $self->status_message("UR_DBI_NO_COMMIT: ".$ENV{UR_DBI_NO_COMMIT});
        UR::Context->commit(); # warning: most code should NEVER do this in a pipeline

        my $instrument_data = Genome::InstrumentData::Imported->get(
            original_data_path => $original_data_path
        );
        unless ($instrument_data) {
            die $self->error_message( "Failed to find new instrument data $original_data_path!");
        }
        if ($instrument_data->__changes__) {
            die "Unsaved changes present on instrument data $instrument_data->{id} from $original_data_path!!!";
        }
        for my $lock (@locks) {
            if ($lock) {
                unless(Genome::Sys->unlock_resource(resource_lock => $lock)) {
                    die $self->error_message("Failed to unlock " . $lock->resource_lock . ".");
                }
            }
        }
        push @instrument_data, $instrument_data;
    }

    return @instrument_data;
}

sub lock {
    my $self = shift;
    my @parts = @_;
    my $lock_key = join('_', @parts);
    $self->debug_message("Creating lock on $lock_key...");
    my $resource_lock = File::Spec->join($ENV{GENOME_LOCK_DIR}, $lock_key);
    my $lock = Genome::Sys->lock_resource(
        resource_lock => $resource_lock,
        max_try => 2,
    );
    return $lock;
};

sub execute_or_die {
    my ($self, $command) = @_;
    my $class = $command->class;

    $self->status_message("Executing $class, command object details are:\n" . Data::Dumper::Dumper($command));

    if ($command->execute){
        $self->status_message("Execution of $class complete.");
        return 1;
    }
    else {
        die $self->error_message("Failed to execute $class.");
    }
}

sub assign_missing_instrument_data_to_model{
    my ($self, $model, @instrument_data) = shift;
}

sub symlink {
    my $self = shift;
    my ($source, $target) = @_;
    if(-l $target && readlink($target) ne $source) {
        die $self->error_message("$target already exists but points to " . readlink($target));
    }
    elsif(! -l $target) {
        Genome::Sys->create_symlink($source, $target);
    }
}

sub get_bam_and_flagstat_from_build{
    my ($self, $build) = @_;
    $self->status_message("getting bam and flagstat from build: ".$build->id);

    #we have to try two different methods to find the merged bam, one uses the <picard software result id>.bam which is the new way, and the old uses the <build_id>.bam

    my @swus = Genome::SoftwareResult::User->get(user_id=>$build->id);

    my @software_results = grep { ref($_) eq "Genome::InstrumentData::AlignmentResult::Merged" } Genome::SoftwareResult->get(id=>[map {$_->software_result_id} @swus]);

    if (scalar @software_results > 1) {
        die $self->error_message("Found more than one merged bam alignment software result for metagenomic alignment build ".$build->id);
    }
    my $alignment = shift @software_results;
    my $alignment_file;
    my $alignment_dir;
    if ($alignment){
        $alignment_dir = $alignment->output_dir;
        $alignment_file = $alignment_dir . "/" . $alignment->id . ".bam";
        unless (-e $alignment_file){
            die $self->error_message("merged bam ($alignment_file) for build ".$build->id." from software result doesn't exist");
        }
    }else{
        #try old method
        $alignment_dir = $build->accumulated_alignments_directory;
        my $build_id = $build->id;
        my @alignment_file = glob("$alignment_dir/$build_id*bam");
        unless (@alignment_file){
            die $self->error_message("no bam file in alignment directory $alignment_dir");
        }
        if (@alignment_file > 1){
            die $self->error_message("more than one bam file found in alignment directory $alignment_dir");
        }
        $alignment_file = shift @alignment_file;
        unless (-e $alignment_file){
            die $self->error_message("Failed to find bam for build ".$build->id." in alignments dir $alignment_dir");
        }
    } 

    my $flagstat_file = "$alignment_file.flagstat";
    unless (-e $flagstat_file){
        die $self->error_message("Failed to flagstat for build ".$build->id." in alignments dir $alignment_dir");
    }
    return ($alignment_file, $flagstat_file);
}

# todo move to the model class
sub build_if_necessary_and_wait {
    my ($self, @models) = @_;

    my (@succeeded_builds, @watched_builds);
    for my $model ( @models ) {
        $self->status_message('Model: '. $model->__display_name__);
        $self->status_message('Search for succeeded build');
        my $succeeded_build = $model->last_succeeded_build;
        if ( $succeeded_build and $self->_verify_model_and_build_instrument_data_match($model, $succeeded_build) ) {
            $self->status_message('Found succeeded build: '.$succeeded_build->__display_name__);
            push @succeeded_builds, $succeeded_build;
            next;
        }
        $self->status_message('No succeeded build');
        $self->status_message('Search for scheduled or running build');
        my $watched_build = $self->_find_scheduled_or_running_build_for_model($model);
        if ( not $watched_build ) {
            $self->status_message('No scheduled or running build');
            $self->status_message('Start build');
            $watched_build = $self->_start_build_for_model($model);
            return if not $watched_build;
        }
        $self->status_message('Watching build: '.$watched_build->__display_name__);
        push @watched_builds, $watched_build;
    }

    my @builds = (@succeeded_builds, @watched_builds);
    if ( not @builds ) {
        $self->error_message('Failed to find or start any builds');
        return;
    }

    if ( @models != @builds ) {
        $self->error_message('Failed to find or start a build for each model');
        return;
    }

    if ( @watched_builds ) {
        for my $build ( @watched_builds ) {
            $self->status_message('Watching build: '.$build->__display_name__);
            $self->wait_for_build($build);
            unless ($build->status eq 'Succeeded') {
                $self->error_message("Failed to execute build (".$build->__display_name__."). Status: ".$build->status);
                return;
            }
        }
    }

    return ( @builds > 1 ? @builds : $builds[0] );
}

sub _verify_model_and_build_instrument_data_match {
    my ($self, $model, $build) = @_;

    Carp::confess('No model to verify instrument data') if not $model;
    Carp::confess('No build to verify instrument data') if not $build;

    my @build_instrument_data = sort {$a->id <=> $b->id} $build->instrument_data;
    my @model_instrument_data = sort {$a->id <=> $b->id} $model->instrument_data;

    $self->status_message('Model: '.$model->__display_name__);
    $self->status_message('Model instrument data: '.join(' ', map { $_->id } @model_instrument_data));
    $self->status_message('Build: '.$build->__display_name__);
    $self->status_message('Build instrument data: '.join(' ', map { $_->id } @build_instrument_data));

    if ( @build_instrument_data != @model_instrument_data ) {
        $self->status_message('Model and build instrument data count does not match');
        return;
    }

    for ( my $i = 0; $i < @model_instrument_data; $i++ ) {
        my $build_instrument_data = $build_instrument_data[$i];
        my $model_instrument_data = $model_instrument_data[$i];

        if ($build_instrument_data->id ne $model_instrument_data->id) {
            $self->status_message("Missing instrument data.");
            return;
        }
    }

    return 1;
}

sub _find_scheduled_or_running_build_for_model {
    my ($self, $model) = @_;

    Carp::confess('No model sent to find running or scheduled build') if not $model;

    $self->status_message('Find running or scheduled build for model: '.$model->__display_name__);

    UR::Context->reload('Genome::Model::Build', model_id => $model->id);

    my $build = $model->latest_build;
    if ( $build and grep { $build->status eq $_ } (qw/ Scheduled Running /) ) {
        return $build;
    }

    return;
}

sub _start_build_for_model {
    my ($self, $model) = @_;

    Carp::confess('No model sent to start build') if not $model;

    my $cmd = 'genome model build start '.$model->id.' --job-dispatch apipe --server-dispatch workflow'; # these are defaults
    $self->status_message('Cmd: '.$cmd);
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv ) {
        $self->error_message('Failed to execute build start command');
        return;
    }

    my $build = $self->_find_scheduled_or_running_build_for_model($model);
    if ( not $build ) {
        $self->error_message('Executed build start command, but cannot find build.');
        return;
    }

    return $build;
}

# TODO: move this to the build class itself
sub wait_for_build {
    my ($self, $build) = @_;
    my $last_status = '';
    my $time = 0;
    my $inc = 30;
    while (1) {
        UR::Context->current->reload($build->the_master_event);
        my $status = $build->status;
        if ($status and !($status eq 'Running' or $status eq 'Scheduled')){
            return 1;
        }

        if ($last_status ne $status or !($time % 300)){
            $self->status_message("Waiting for build(~$time sec) ".$build->id.", status: $status");
        }
        sleep $inc;
        $time += 30;
        $last_status = $status;
    }
}

# TODO: move this to the model class itself
#this name is maybe not the best
sub _process_sra_instrument_data {
    my ($self, $instrument_data) = @_;
    my $lane = $instrument_data->lane;

    my %params;
    $params{instrument_data} = $instrument_data;
    $params{n_removal_threshold} = $self->n_removal_threshold if $self->n_removal_threshold;
    $params{non_n_base_threshold} = $self->non_n_base_threshold if $self->non_n_base_threshold;
    $params{dust} = $self->dust_unaligned_reads if $self->dust_unaligned_reads;

    my $cmd = Genome::InstrumentData::Command::PostProcessAndImport->create(%params);
    unless($cmd){
        die $self->error_message("Couldn't create PostProcessAndImport command for instrument data ".$instrument_data->id);
    }
    my $rv = $cmd->execute;
    unless ($rv){
        die $self->error_message("Couldn't execute PostProcessAndImport command for instrument data ".$instrument_data->id);
    }

    my @instrument_data = $cmd->post_processed_instrument_data();

    unless (@instrument_data){
        die $self->error_message("No post-processed instrument data returned as output from PostProcessAndImport command for instrument data ". $instrument_data->id);
    }

    return @instrument_data;
}

sub _process_unaligned_reads {
    my ($self, $alignment) = @_;

    my $instrument_data = $alignment->instrument_data;
    my $lane = $instrument_data->lane;
    my $instrument_data_id = $instrument_data->id;

    my $dir = $alignment->output_dir;
    my $bam = $dir . '/all_sequences.bam';
    unless (-e $bam) {
        $self->error_message("Failed to find expected BAM file $bam\n");
        return;
    }

    my $tmp_dir = "$UNALIGNED_TEMPDIR/unaligned_reads";
    $tmp_dir .= "/".$alignment->id;
    my @instrument_data;
    {
        my $subdir = 'n-remove_'.$self->n_removal_threshold; # FIXME uninit warnings
        if ($self->dust_unaligned_reads){
            $subdir.='/dusted';
        }
        $self->status_message("Preparing imported instrument data for import path $tmp_dir/$subdir");

        # proceed extracting and uploading unaligned reads into $tmp_dir/$subdir....
        # resolve the paths at which we will place processed instrument data
        # we're currently using these paths to find previous unaligned reads processed the same way

        my $forward_basename = "s_$lane" . "_1_sequence.txt";
        my $reverse_basename = "s_$lane" . "_2_sequence.txt";
        my $fragment_basename = "s_$lane" . "_sequence.txt";

        my $forward_unaligned_data_path     = "$tmp_dir/$instrument_data_id/$forward_basename";
        my $reverse_unaligned_data_path     = "$tmp_dir/$instrument_data_id/$reverse_basename";
        my $fragment_unaligned_data_path    = "$tmp_dir/$instrument_data_id/$fragment_basename";

        my @expected_original_paths;
        my $expected_data_path0 = "$tmp_dir/$subdir/$fragment_basename";
        my $expected_data_path1 = "$tmp_dir/$subdir/$forward_basename";
        my $expected_data_path2 = "$tmp_dir/$subdir/$reverse_basename";


        my $expected_se_path = $expected_data_path0;
        my $expected_pe_path = $expected_data_path1 . ',' . $expected_data_path2;


        my @upload_paths;
        my ($se_lock, $pe_lock);


        # check for previous unaligned reads
        $self->status_message("Checking for previously imported unaligned and post-processed reads from: $tmp_dir/$subdir");
        my $se_instdata = Genome::InstrumentData::Imported->get(original_data_path => $expected_se_path);
        if ($se_instdata) {
            $self->status_message("imported instrument data already found for path $expected_se_path, skipping");
        }
        else {
            $se_lock = $self->lock($instrument_data_id, basename($expected_se_path));
            unless ($se_lock) {
                die $self->error_message("Failed to lock $expected_se_path.");
            }
            push @upload_paths, $expected_se_path;
        }

        my $pe_instdata = Genome::InstrumentData::Imported->get(original_data_path => $expected_pe_path);
        if ($pe_instdata) {
            $self->status_message("imported instrument data already found for path $expected_pe_path, skipping");
        }
        elsif ( $instrument_data->is_paired_end ) {
            $pe_lock = $self->lock($instrument_data_id, basename($expected_pe_path));
            unless ($pe_lock) {
                die $self->error_message("Failed to lock $expected_pe_path.");
            }
            push @upload_paths, $expected_pe_path;
        }

        unless (@upload_paths) {
            $self->status_message("skipping read processing since all data is already processed and uploaded");
            return grep { defined } ($se_instdata, $pe_instdata);
        }

        for my $path (@upload_paths) {
            $self->status_message("planning to upload for $path");
        }

        # extract
        my $create_directory = eval{ Genome::Sys->create_directory($tmp_dir.'/'.$subdir); };
        if ( not $create_directory ) {
            die "Failed to create tmp directory ($tmp_dir/$subdir): $@";
        }
        $self->status_message("Preparing imported instrument data for import path $tmp_dir/$subdir");
        $self->status_message("Extracting unaligned reads from $bam");
        my $cmd = "genome-perl5.10 -S gmt bio-samtools bam-to-unaligned-fastq --bam-file $bam --output-directory $tmp_dir";
        my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
        if ( not $rv ) {
            die "Failed to extract unaligned reads: $@";
        }

        #$self->status_message("Perl version: $]");
        #my $extract_unaligned = Genome::Model::Tools::BioSamtools::BamToUnalignedFastq->create(
        #    bam_file => $bam,
        #    output_directory =>$tmp_dir,
        #);
        #my $rv = $extract_unaligned->execute;
        #unless ($rv){
        #    die $self->error_message("Couldn't extract unaligned reads from bam file $bam");
        #}

        my @expected_output_fastqs = ( $instrument_data->is_paired_end )
        ?  ($forward_unaligned_data_path, $reverse_unaligned_data_path, $fragment_unaligned_data_path)
        :  ($fragment_unaligned_data_path);
        my @missing = grep {! -e $_} grep { defined($_) and length($_) } @expected_output_fastqs;
        if (@missing){
            die $self->error_message(join(", ", @missing)." unaligned files missing after bam extraction");
        }
        $self->status_message("Extracted unaligned reads from bam file (@expected_output_fastqs)");

        # process the fragment data
        if (-e $fragment_unaligned_data_path) {
            $self->status_message("processind single-end reads...");
            my $processed_fastq = $self->_process_unaligned_fastq($fragment_unaligned_data_path, $expected_data_path0);
            unless (-e $expected_data_path0){
                die $self->error_message("Expected data path does not exist after fastq processing: $expected_data_path0");
            }
        }

        # process the paired data
        if (-e $forward_unaligned_data_path or -e $reverse_unaligned_data_path) {
            $self->status_message("processind paired-end reads...");
            unless (-e $forward_unaligned_data_path and -e $reverse_unaligned_data_path) {
                die "Missing forward and reverse unaligned data?";
            }
            my $processed_fastq1 = $self->_process_unaligned_fastq($forward_unaligned_data_path, $expected_data_path1);
            my $processed_fastq2 = $self->_process_unaligned_fastq($reverse_unaligned_data_path, $expected_data_path2);
            my @missing = grep {! -e $_} ($expected_data_path1, $expected_data_path2);
            if (@missing){
                $self->error_message("Expected data paths do not exist after fastq processing: ".join(", ", @missing));
                Carp::confess($self->error_message);
            }
        }

        # upload
        $self->status_message("uploading new instrument data from the post-processed unaligned reads...");
        my @properties_from_prior = qw/
        run_name
        sequencing_platform
        median_insert_size
        sd_above_insert_size
        library_name
        sample_name
        /;
        my @errors;
        my %properties_from_prior;
        for my $property_name (@properties_from_prior) {
            my $value = $instrument_data->$property_name;
            no warnings;
            $self->status_message("Value for $property_name is $value");
            $properties_from_prior{$property_name} = $value;
        }
        $properties_from_prior{subset_name} = $instrument_data->lane;

        #my @instrument_data;
        for my $original_data_path (@upload_paths) {
            $self->status_message("Attempting to upload $original_data_path...");
            if ($original_data_path =~ /,/){
                $properties_from_prior{is_paired_end} = 1;
            }else{
                $properties_from_prior{is_paired_end} = 0;
            }
            #my $previous = Genome::InstrumentData::Imported->get(
            #    original_data_path => $original_data_path,
            #);
            my $previous;
            if ($previous){
                $self->error_message("imported instrument data already found for path $original_data_path????");
                Carp::confess($self->error_message);
                #push @instrument_data, $previous;
                #next;
            }
            my %params = (
                %properties_from_prior,
                source_data_files => $original_data_path,
                import_format => 'illumina fastq',
            );
            $self->status_message("importing fastq with the following params:" . Data::Dumper::Dumper(\%params));


            my $command = Genome::InstrumentData::Command::Import::Fastq->create(%params);
            unless ($command) {
                $self->error_message( "Couldn't create command to import unaligned fastq instrument data!");
            };
            my $result = $command->execute();
            unless ($result) {
                die $self->error_message( "Error importing data from $original_data_path! " . Genome::InstrumentData::Command::Import::Fastq->error_message() );
            }
            $self->status_message("committing newly created imported instrument data");
            $self->status_message("UR_DBI_NO_COMMIT: ".$ENV{UR_DBI_NO_COMMIT});
            UR::Context->commit(); # warning: most code should NEVER do this in a pipeline

            my $new_instrument_data = Genome::InstrumentData::Imported->get($command->generated_instrument_data_id);
            unless ($new_instrument_data) {
                die $self->error_message( "Failed to find new instrument data $original_data_path!");
            }
            if ($new_instrument_data->__changes__) {
                die "Unsaved changes present on instrument data $new_instrument_data->{id} from $original_data_path!!!";
            }
            if (!$new_instrument_data->is_paired_end && $se_lock) {
                $self->status_message("Attempting to remove lock on $se_lock...");
                unless(Genome::Sys->unlock_resource(resource_lock => $se_lock)) {
                    die $self->error_message("Failed to unlock $se_lock.");
                }
                undef($se_lock);
            }
            if ($new_instrument_data->is_paired_end && $pe_lock) {
                $self->status_message("Attempting to remove lock on $pe_lock...");
                unless(Genome::Sys->unlock_resource(resource_lock => $pe_lock)) {
                    die $self->error_message("Failed to unlock $pe_lock.");
                }
                undef($pe_lock);
            }
            push @instrument_data, $new_instrument_data;
        }
        #return @instrument_data;
    };

    system "/bin/rm -rf $tmp_dir";

    if ( not @instrument_data ) {
        die $self->error_message("Error processing unaligned reads!");
    }

    return @instrument_data;
}

sub _process_unaligned_fastq_pair {
    my $self = shift;
    my ($forward, $reverse, $forward_out, $reverse_out, $fragment_out) = @_;
    #run dust on forward and reverse
    my $forward_dusted;
    my $reverse_dusted;

    if ($self->dust_unaligned_reads){
        $self->status_message("Dusting fastq pair $forward, $reverse");
        $forward_dusted = "$forward.DUSTED";
        $reverse_dusted = "$reverse.DUSTED";

        $self->dust_fastq($forward, $forward_dusted);
        $self->dust_fastq($reverse, $reverse_dusted);
    }else{
        $self->status_message("skipping dusting");
        $forward_dusted = $forward;
        $reverse_dusted = $reverse;
    }

    #run pairwise n-removal
    if ($self->n_removal_threshold){
        $self->status_message("running remove-n-pairwise on $forward, $reverse");
        my $cmd = Genome::Model::Tools::Fastq::RemoveNPairwise->create(
            forward_fastq => $forward_dusted,
            reverse_fastq => $reverse_dusted,
            forward_n_removed_file => $forward_out,
            reverse_n_removed_file => $reverse_out,
            singleton_n_removed_file => $fragment_out,
            n_removal_threshold => $self->n_removal_threshold,
        );
        unless ($cmd){
            die $self->error_message("couldn't create remove-n-pairwise command for $forward_dusted, $reverse_dusted!");
        }
        my $rv = $cmd->execute;
        unless ($rv){
            die $self->error_message("couldn't create remove-n-pairwise command for $forward_dusted, $reverse_dusted!");
        }
        unless(-e $forward_out && -e $reverse_out && -e $fragment_out){
            die $self->error_message("couldn't find all expected output files! $forward_out, $reverse_out, $fragment_out");
        }
        #clean up, maybe make these temp files
        if ($self->dust_unaligned_reads){
            #only need to do this if we actually dusted
            unlink $forward_dusted;
            unlink $reverse_dusted;
        }

        #return the 3 processed fastq files
        return ($forward_out, $reverse_out, $fragment_out);
    }else{
        $self->status_message("skipping n-removal");
        Genome::Sys::copy_file($forward_dusted, $forward_out);
        Genome::Sys::copy_file($reverse_dusted, $reverse_out);
        if ($self->dust_unaligned_reads){
            #only need to do this if we actually dusted
            unlink $forward_dusted;
            unlink $reverse_dusted;
        }
        return ($forward_out, $reverse_out);
    }
}

sub dust_fastq{
    my ($self, $in, $out) = @_;
    my $cmd = Genome::Model::Tools::Fastq::Dust->create(
        fastq_file => $in,
        output_file => $out,
    );
    unless ($cmd){
        die $self->error_message("couldn't create dust command for $in -> $out!");
    }
    my $rv = $cmd->execute;
    unless ($rv){
        die $self->error_message("failed to execute dust command for $in -> $out! rv:$rv");
    }
    unless (-s $out){
        die $self->error_message("expected output file $out doesn't exist or has 0 size!");
    }
    return $out;
}

sub _process_unaligned_fastq {
    my $self = shift;
    my ($fastq_file, $output_path) = @_;

    my $dusted_fastq;
    if ($self->dust_unaligned_reads){
        $dusted_fastq = "$fastq_file.DUSTED";
        $self->dust_fastq($fastq_file, $dusted_fastq);
    }else{
        $self->status_message("skipping dusting $fastq_file");
        $dusted_fastq = $fastq_file;
    }

    if ($self->n_removal_threshold){
        $self->status_message("Running n-removal on file $fastq_file");
        my $cmd = Genome::Model::Tools::Fastq::RemoveN->create(
            fastq_file => $dusted_fastq,
            n_removed_file => $output_path,
            n_removal_threshold => $self->n_removal_threshold,
        );
        unless ($cmd){
            die $self->error_message("couldn't create remove-n command for $dusted_fastq");
        }
        my $rv = $cmd->execute;
        unless ($rv){
            die $self->error_message("couldn't execute remove-n command for $dusted_fastq");
        }
    } else {
        $self->status_message("No n-removal cutoff specified, skipping");
        Genome::Sys->copy_file($dusted_fastq, $output_path);
    }
    if ($self->dust_unaligned_reads){
        unlink $dusted_fastq;
    }
    return $output_path;
}




1;
