
package Genome::Model::ReferenceAlignment;

#:eclark 11/18/2009 Code review.

# I'm not sure that we need to have all these via properties.  Does it really gain us that much code clarity?
# Some of the other properties seem generic enough to be part of Genome::Model, not this subclass.
# The entire set of _calculcute* methods could be refactored away.
# Several deprecated/todo comments scattered in the code below that should either be removed or implemented.
# Most of the methods at the bottom are for reading/writing of a gold_snp_file, this should be implemented as
# part of Genome/Utility/IO

use strict;
use warnings;

use Genome;
use Term::ANSIColor;
use File::Path;
use File::Basename;
use File::Spec;
use IO::File;
use Sort::Naturally;
use Data::Dumper;

my %DEPENDENT_PROPERTIES = (
    # dbsnp_build and annotation depend on reference_sequence_build
    'reference_sequence_build' => [
        'dbsnp_build',
        'annotation_reference_build',
    ],
);

class Genome::Model::ReferenceAlignment {
    is => 'Genome::ModelDeprecated',
    has_mutable_input => [
        genotype_microarray => {
            is => 'Genome::Model::GenotypeMicroarray',
            is_optional => 1,
            doc => 'Genotype Microarray model used for QC and Gold SNP Concordance report', 
            # this is redundant with genotype_microarray_model, which has the correct
            # method name, but uses "genotype_microarray" in the db layer (fix that)
        },
        reference_sequence_build => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            doc => 'reference sequence to align against',
        },
        dbsnp_build => {
            is => 'Genome::Model::Build::ImportedVariationList',
            doc => 'dbsnp build to compare against',
            is_optional => 1,
        },
        annotation_reference_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
            doc => 'The reference build used for variant annotation',
            is_optional => 1,
        },
    ],
    has_param => [
        sequencing_platform => {
            doc => 'The sequencing platform from whence the model data was generated',
            valid_values => ['454', 'solexa', 'sanger'],
        },
        dna_type => {
            doc => 'the type of dna used in the reads for this model',
            valid_values => ['genomic dna', 'cdna']
        },
        transcript_variant_annotator_version => {
            doc => 'Version of the "annotate transcript-variants" tool to run during the annotation step',
            is_optional => 1,
            default_value => Genome::Model::Tools::Annotate::TranscriptVariants->default_annotator_version,
            valid_values => [ 0,1,2,3,4],#Genome::Model::Tools::Annotate::TranscriptVariants->available_versions ],
        },
        transcript_variant_annotator_filter => {
            doc => 'annotation-filter option to be used by the "annotate transcript-variants" tool run during the annotation step',
            is_optional => 1,
            default_value => 'top',
            valid_values => ['top', 'none', 'gene'],
        },
        transcript_variant_annotator_accept_reference_IUB_codes => {
            doc => 'annotation accept-reference-IUB-codes option to be used by the "annotate transcript-variants" to run during the annotation step',
            is_optional => 1,
            default_value => '0',
            valid_values => [0, 1],
        },
        alignment_strategy => {
            is => 'Text',
            is_many => 0,
            is_optional => 1,
            doc => 'Strategy to be used to align the instrument data (this will eventually replace many of the alignment parameters)',
        },
        snv_detection_strategy => {
            is => "Text",
            is_many => 0,
            is_optional =>1,
            doc => "Strategy to be used to detect snvs.",
        },
        indel_detection_strategy => {
            is => "Text",
            is_many => 0,
            is_optional =>1,
            doc => "Strategy to be used to detect indels.",
        },
        sv_detection_strategy => {
            is => "Text",
            is_many => 0,
            is_optional =>1,
            doc => "Strategy to be used to detect svs.",
        },
        cnv_detection_strategy => {
            is => "Text",
            is_many => 0,
            is_optional =>1,
            doc => "Strategy to be used to detect cnvs.",
        },
        picard_version => {
            doc => 'picard version for MarkDuplicates, MergeSamfiles, CreateSequenceDictionary...',
            is_optional => 1,
        },
        samtools_version => {
            doc => 'samtools version for SamToBam, samtools merge, etc...',
            is_optional => 1,
        },
        bedtools_version => {
            doc => 'bedtools version for bedtools bamtofastq',
            is_optional => 1,
        },
        merger_name => {
            doc => 'name of bam merger, picard, samtools (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        merger_version => {
            doc => 'version of bam merger (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        merger_params => {
            doc => 'parameters of bam merger (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        duplication_handler_name => {
            doc => 'name of pcr duplication handler (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        duplication_handler_version => {
            doc => 'version of pcr duplication handler (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        duplication_handler_params => {
            doc => 'parameters of pcr duplication handler (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        read_aligner_name => {
            doc => 'alignment algorithm/software used for this model (this will be replaced by alignment_strategy)',
        },
        read_aligner_version => {
            doc => 'the aligner version used for this model (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        read_aligner_params => {
            doc => 'command line args for the aligner (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        force_fragment => {
            is => 'Integer',
            #This doesn't seem to work yet because of the create code, can't the valid values logic be removed from create???
            default_value => '0',
            #valid_values => ['0', '1'],
            doc => 'force all alignments as fragment reads',
            is_optional => 1,
        },
        read_trimmer_name => {
            doc => 'trimmer algorithm/software used for this model',
            is_optional => 1,
        },
        read_trimmer_version => {
            doc => 'the trimmer version used for this model',
            is_optional => 1,
        },
        read_trimmer_params => {
            doc => 'command line args for the trimmer',
            is_optional => 1,
        },
        coverage_stats_params => {
            doc => 'parameters necessary for generating reference coverage in the form of two comma delimited lists split by a colon like 1,5,10,15,20:0,200,500',
            is_optional => 1,
        },
        append_event_steps => {
            doc => 'Event classes to append to event_stage_job_classes, e.g. "alignment => Genome::Model::Event::Build::ReferenceAlignment::QC::CopyNumber".',
            is_optional => 1,
            is_deprecated => 1,
        },
    ],
    has => [
        subject => {
            is => 'Genome::Sample',
            id_by => 'subject_id',
            doc => 'the subject of alignment and variant detection is a single sample'
        },
        dbsnp_model => {
            via => 'dbsnp_build',
            to => 'model',
            is_mutable => 1,
        },

        # TODO: these are the right accessors, but the underlying input name is wrong :(
        # fix the db, rename the above genotype_microarray to have _model, and get rid of these
        genotype_microarray_model_id => {
            is => 'Text',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'genotype_microarray', 'value_class_name' => 'Genome::Model::GenotypeMicroarray', ],
            is_mutable => 1,
            is_optional => 1,
            doc => 'Genotype Microarray model used for QC and Gold SNP Concordance report',
        },
        genotype_microarray_model => {
            is => 'Genome::Model::GenotypeMicroarray',
            id_by => 'genotype_microarray_model_id',
        },

        reference_sequence_name      => { via => 'reference_sequence_build', to => 'name' },
        annotation_reference_name    => { via => 'annotation_reference_build', to => 'name' },
        target_region_set_name => {
            is_many => 1, is_mutable => 1, is => 'Text', via => 'inputs', to => 'value_id', 
            where => [ name => 'target_region_set_name', value_class_name => 'UR::Value' ],
            is_optional => 1,
        },
    ],
    doc => 'A genome model produced by aligning DNA reads to a reference sequence.'
};

sub create {
    my $class = shift;

    # This is a temporary hack to allow annotation_reference_build (currently calculated) to be
    # passed in as an object. Once the transition to using model inputs for this parameter vs
    # processing profile params, annotation_reference_build can work like reference_sequence_build
    # and this code can go away.
    my @args = @_;
    if (scalar(@_) % 2 == 0) {
        my %args = @args;
        if (defined $args{annotation_reference_build}) {
            $args{annotation_reference_build_id} = (delete $args{annotation_reference_build})->id;
            @args = %args;
        }
    }

    my $self = $class->SUPER::create(@args)
        or return;

    unless ( $self->reference_sequence_build ) {
        $self->error_message("Missing needed reference sequence build during reference alignment model creation.");
        $self->delete;
        return;
    }

    unless ($self->genotype_microarray_model_id) {
        my $genotype_model = $self->default_genotype_model;
        if ($genotype_model) {
            $self->genotype_microarray_model_id($genotype_model->id);
        }
    }

    if ($self->read_aligner_name and $self->read_aligner_name eq 'newbler') {
        my $new_mapping = Genome::Model::Tools::454::Newbler::NewMapping->create(
            dir => $self->alignments_directory,
        );
        unless ($self->new_mapping) {
            $self->error_message('Could not setup newMapping for newbler in directory '. $self->alignments_directory);
            return;
        }
        my @fasta_files = grep {$_ !~ /all_sequences/} $self->get_subreference_paths(reference_extension => 'fasta');
        my $set_ref = Genome::Model::Tools::454::Newbler::SetRef->create(
                                                                    dir => $self->alignments_directory,
                                                                    reference_fasta_files => \@fasta_files,
                                                                );
        unless ($set_ref->execute) {
            $self->error_message('Could not set refrence setRef for newbler in directory '. $self->alignments_directory);
            return;
        }
    }
    return $self;
}

sub region_of_interest_set {
    my $self = shift;

    my $name = $self->region_of_interest_set_name;
    return unless $name;
    my $roi_set = Genome::FeatureList->get(name => $name);
    unless ($roi_set) {
        die('Failed to find feature-list with name: '. $name);
    }
    return $roi_set;
}

sub accumulated_alignments_directory {
    my $self = shift;
    my $last_complete_build = $self->last_complete_build;
    return if not $last_complete_build;
    return File::Spec->join($last_complete_build->data_directory, 'alignments');
}

sub is_capture {
    my $self = shift;
    if (defined $self->target_region_set_name) {
        return 1;
    }
    return 0;
}

sub is_lane_qc {
    my $self = shift;
    my $pp = $self->processing_profile;
    if ($pp->append_event_steps && $pp->append_event_steps =~ /LaneQc/) {
        return 1;
    }
    return 0;
}

sub build_needed {
    my $self = shift;
    my $needed = $self->SUPER::build_needed;

    if($self->is_lane_qc) {
        $needed &&= $self->genotype_microarray and $self->genotype_microarray->last_complete_build;
    }

    return $needed;
}

# Determines the correct genotype model to use via the official genotype data assigned to the subject
sub default_genotype_model {
    my $self = shift;
    my $sample = $self->subject;
    return unless $sample->isa('Genome::Sample');

    my @genotype_models = sort { $a->creation_date cmp $b->creation_date } $sample->default_genotype_models;
    return unless @genotype_models;

    @genotype_models = grep { $_->reference_sequence_build && $_->reference_sequence_build->is_compatible_with($self->reference_sequence_build) } @genotype_models;
    return unless @genotype_models;


    my $chosen_genotype_model;
    if (@genotype_models > 1) {
        my $message = "Found multiple compatible genotype models for sample " . $sample->id . " and reference alignment model " . $self->id;
        my @used_genotype_models = grep { $_->is_used_as_model_or_build_input or $_->builds_are_used_as_model_or_build_input } @genotype_models;
        if (@used_genotype_models) {
            $self->warning_message($message . ", choosing the most recent previously-used model.");
            $chosen_genotype_model = $used_genotype_models[-1];
        } else {
            $self->warning_message($message . ", choosing most recent model.");
            $chosen_genotype_model = $genotype_models[-1];
        }
    }
    else {
        $chosen_genotype_model = $genotype_models[0];
    }

    return $chosen_genotype_model;
}

sub dependent_properties {
    my ($self, $property_name) = @_;
    return @{$DEPENDENT_PROPERTIES{$property_name}} if exists $DEPENDENT_PROPERTIES{$property_name};
    return;
}

sub check_for_updates {
    # TODO: make an observer in ::Site::TGI and move the method below and its kin there.
    # It should watch the "create" signal for Genome::Model::Build::ReferenceAlignment.
    my $self = shift;
    $self->check_and_update_genotype_input;
    return 1;
}

sub check_and_update_genotype_input {
    my $self = shift;
    my $default_genotype_model = $self->default_genotype_model;
    return 1 unless $default_genotype_model;

    if (defined $self->genotype_microarray_model_id and $self->genotype_microarray_model_id ne $default_genotype_model->id) {
        if(Genome::Sys->user_has_role($self->run_as, 'production')) {
            $self->warning_message("Sample " . $self->subject_id . " points to genotype model " . $default_genotype_model->id .
                ", which disagrees with the genotype model input of this model (" . $self->genotype_microarray_model_id .
                "), overwriting!");
            $self->genotype_microarray_model_id($default_genotype_model->id);
        }
    }
    elsif (not defined $self->genotype_microarray_model_id) {
        $self->genotype_microarray_model_id($default_genotype_model->id);
    }

    return 1;
}


sub default_lane_qc_model_name_for_instrument_data {
    my $self = shift;
    my $instrument_data = shift;

    my $subset_name = $instrument_data->subset_name || 'unknown-subset';
    my $run_name_method = $instrument_data->can('short_name') ? 'short_name' : 'run_name';
    my $run_name = $instrument_data->$run_name_method || 'unknown-run';
    my $lane_name = $run_name . '.' . $subset_name;
    my $model_name = $lane_name . '.prod-qc';

    if ($instrument_data->target_region_set_name) {
        $model_name .= '.capture.' . $instrument_data->target_region_set_name;
    }

    return $model_name;
}


sub default_lane_qc_model_for_instrument_data {
    my $class           = shift;
    my @instrument_data = @_;

    my @lane_qc_models;
    for my $instrument_data (@instrument_data) {
        my $lane_qc_model_name = $class->default_lane_qc_model_name_for_instrument_data($instrument_data);
        my $lane_qc_model = Genome::Model::ReferenceAlignment->get(name => $lane_qc_model_name);
        push @lane_qc_models, $lane_qc_model if $lane_qc_model;
    }

    return @lane_qc_models;
}

sub qc_type_for_target_region_set_name {
    my $class = shift;
    my $target_region_set_name = shift;
    return ($target_region_set_name ? 'capture' : 'wgs');
}

# FIXME This needs to be renamed/refactored. The method name does not accurately describe what
# this method actually does.
sub get_lane_qc_models {
    my $self = shift;

    my $subject = $self->subject;

    unless ($subject->default_genotype_data_id) {
        $self->warning_message("Sample is missing default_genotype_data_id cannot create lane QC model.");
        return;
    }

    my @lane_qc_models;
    my @instrument_data = sort { $a->run_name . $a->subset_name cmp $b->run_name . $b->subset_name } $self->instrument_data;
    for my $instrument_data (@instrument_data) {
        my $lane_qc_model_name = $self->default_lane_qc_model_name_for_instrument_data($instrument_data);

        my $existing_model = Genome::Model->get(name => $lane_qc_model_name);
        my $has_current_qc =
            $existing_model &&
            $existing_model->config_profile_item &&
            $existing_model->config_profile_item->is_current;
        if ($has_current_qc) {
            $self->debug_message("Default lane QC model " . $existing_model->__display_name__ . " already exists.");

            my @existing_instrument_data = $existing_model->instrument_data;
            unless (@existing_instrument_data) {
                $existing_model->add_instrument_data($instrument_data);
                $existing_model->build_requested(1, 'instrument data assigned');
                $self->debug_message("New build requested for lane qc model " . $existing_model->__display_name__ .
                    " because it just had instrument data assigned to it"
                );
            }

            unless ($existing_model->genotype_microarray_model_id) {
                $self->build_requested(1, 'genotype ' . $subject->default_genotype_data_id . ' data added');
                $self->debug_message("New build requested for lane QC model " . $existing_model->__display_name__ . 
                    " because it is missing the genotype_microarray input."
                );
            }

            push @lane_qc_models, $existing_model;
        }
    }

    return @lane_qc_models;
}

sub _additional_parts_for_default_name {
    my $self = shift;
    my @parts;

    my @regions = $self->target_region_set_name;
    push @parts, 'capture' if @regions;
    push @parts, $self->region_of_interest_set_name if $self->region_of_interest_set_name;

    return @parts;
}

sub default_model_name {
    my $self = shift;

    if ($self->is_lane_qc) {
        die $self->error_message('Attempting to get the default name when creating a lane qc model that does not have instrument data assigned. Please check the analysis project configuration or include the instrument data when creating the model.') if not $self->instrument_data;
        return $self->_get_incremented_name($self->default_lane_qc_model_name_for_instrument_data($self->instrument_data));
    } else {
        return $self->SUPER::default_model_name();
    }
}

sub experimental_subject {
    my $self = shift;

    my $subject = $self->subject;

    return $subject if $subject->isa('Genome::Sample');
    return;
}

sub control_subject {
    my $self = shift;

    return; #reference alignment has no control subject
}

sub _initialize_build {
    my($self,$build) = @_;

    # Check that the annotator version param is sane before doing the build
    my $annotator_version;
    my $worked = eval {
        my $model = $build->model;
        my $pp = $model->processing_profile;
        $annotator_version = $pp->transcript_variant_annotator_version;
        # When all processing profiles have a param for this, remove this unless block so
        # they'll fail if it's missing
        unless (defined $annotator_version) {
            $annotator_version = Genome::Model::Tools::Annotate::TranscriptVariants->default_annotator_version;
        }

        my %available_versions = map { $_ => 1 } Genome::Model::Tools::Annotate::TranscriptVariants->available_versions;
        unless ($available_versions{$annotator_version}) {
            die "Requested annotator version ($annotator_version) is not in the list of available versions: "
                . join(', ',keys(%available_versions));
        }
        1;
    };
    unless ($worked) {
        $self->error_message("Could not determine which version of the Transcript Variants annotator to use: $@");
        return;
    }
    return 1;
}

sub _resolve_workflow_for_build {
    my $self = shift;
    my $build = shift;

    if ($build->processing_profile->append_event_steps) {
        $self->fatal_message('Appended events are no longer supported in this pipeline.');
    }

    my $workflow = Genome::WorkflowBuilder::DAG->create(
        name => $build->workflow_name,
    );

    my $alignment = Genome::WorkflowBuilder::Command->create(
        name => 'alignment',
        command => 'Genome::Model::ReferenceAlignment::Command::AlignReads',
    );
    $workflow->add_operation($alignment);
    $workflow->connect_input(
        input_property => 'build',
        destination => $alignment,
        destination_property => 'build',
    );

    my $coverage = Genome::WorkflowBuilder::Command->create(
        name => 'coverage',
        command => 'Genome::Model::ReferenceAlignment::Command::CoverageStats',
    );
    $workflow->add_operation($coverage);

    my $bam_qc = Genome::WorkflowBuilder::Command->create(
        name => 'bam-qc',
        command => 'Genome::Model::ReferenceAlignment::Command::BamQc',
        parallel_by => 'alignment_result',
    );
    $workflow->add_operation($bam_qc);
    $workflow->create_link(
        source => $alignment,
        source_property => 'individual_alignment_results',
        destination => $bam_qc,
        destination_property => 'alignment_result',
    );

    my $detect_variants = Genome::WorkflowBuilder::Command->create(
        name => 'detect-variants',
        command => 'Genome::Model::ReferenceAlignment::Command::DetectVariants',
    );
    $workflow->add_operation($detect_variants);

    for my $post_alignment_step ($coverage, $bam_qc, $detect_variants) {
        $workflow->create_link(
            source => $alignment,
            source_property => 'build',
            destination => $post_alignment_step,
            destination_property => 'build',
        );
    }

    my $annotation = Genome::WorkflowBuilder::Command->create(
        name => 'annotation',
        command => 'Genome::Model::ReferenceAlignment::Command::AnnotateVariants',
    );
    $workflow->add_operation($annotation);
    $workflow->create_link(
        source => $detect_variants,
        source_property => 'build',
        destination => $annotation,
        destination_property => 'build',
    );

    my @block_inputs = ($annotation, $bam_qc, $coverage);
    my $block = Genome::WorkflowBuilder::Block->create(
        name => 'prepare-for-reports',
        properties => ['build', (map { $_->name } @block_inputs)],
    );
    $workflow->add_operation($block);
    $workflow->connect_input(
        input_property => 'build',
        destination => $block,
        destination_property => 'build',
    );
    for my $block_input (@block_inputs) {
        $workflow->create_link(
            source => $block_input,
            source_property => 'result',
            destination => $block,
            destination_property => $block_input->name,
        );
    }

    my $reports = Genome::WorkflowBuilder::Command->create(
        name => 'reports',
        command => 'Genome::Model::ReferenceAlignment::Command::RunReports',
    );
    $workflow->add_operation($reports);
    $workflow->create_link(
        source => $block,
        source_property => 'build',
        destination => $reports,
        destination_property => 'build',
    );
    $workflow->connect_output(
        source => $reports,
        source_property => 'result',
        output_property => 'report_result',
    );

    $workflow->recursively_set_log_dir($build->log_directory);
    return $workflow;
}

sub map_workflow_inputs {
    my $self = shift;
    my $build = shift;

    return (build => $build);
}

sub individual {
    my $self = shift;
    return $self->subject->individual;
}

1;
