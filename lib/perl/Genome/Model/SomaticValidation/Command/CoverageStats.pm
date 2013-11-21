package Genome::Model::SomaticValidation::Command::CoverageStats;

use strict;
use warnings;

use File::Path qw(rmtree);

use Genome;

class Genome::Model::SomaticValidation::Command::CoverageStats {
    is => ['Command::V2'],
    has_input => [
        build_id => {
            is => 'Number',
            doc => 'ID of the SomaticValidation build upon which to run coverage stats',
        },
        mode => {
            is => 'Text',
            valid_values => ['tumor', 'normal'],
        },
    ],
    has => [
        build => {
            is => 'Genome::Model::Build::SomaticValidation',
            id_by => 'build_id',
            doc => 'The build upon which to run coverage stats',
        },
    ],
    has_param => [
        lsf_queue => {
            default => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
        },
    ],
    has_optional_output => [
        reference_coverage_result_id => {
            is => 'Number',
            doc => 'ID of the result from running this step',
        },
    ],
    has_optional => [
        reference_coverage_result => {
            is => 'Genome::InstrumentData::AlignmentResult::Merged::CoverageStats',
            id_by => 'reference_coverage_result_id',
            doc => 'The result from running this step',
        },
    ],
    doc => 'runs ref-cov on the bam(s) produced in the alignment step',
};

sub sub_command_category { 'pipeline steps' }

sub shortcut {
    my $self = shift;

    unless($self->should_run) {
        $self->status_message('Sample not specified on build; skipping.');
        return 1;
    }

    my %params = $self->params_for_result;
    my $result = Genome::InstrumentData::AlignmentResult::Merged::CoverageStats->get_with_lock(%params);

    if($result) {
        $self->status_message('Using existing result ' . $result->__display_name__);
        return $self->link_result_to_build($result);
    } else {
        return;
    }
}

sub execute {
    my $self = shift;
    my $build = $self->build;

    unless($self->should_run) {
        $self->status_message('Sample not specified on build; skipping.');
        return 1;
    }

    unless($self->_reference_sequence_matches) {
        die $self->error_message;
    }

    my %params = (
        $self->params_for_result,
        log_directory => $build->log_directory,
    );

    my $result = Genome::InstrumentData::AlignmentResult::Merged::CoverageStats->get_or_create(%params);

    $self->link_result_to_build($result);

    my $as_ref = $result->alignment_summary_hash_ref;
    unless ($as_ref) {
        $self->error_message('Failed to load the alignment summary metrics!');
        die($self->error_message);
    }
    my $cov_ref = $result->coverage_stats_summary_hash_ref;
    unless ($cov_ref) {
        $self->error_message('Failed to load the coverage summary metrics!');
        die($self->error_message);
    }

    return 1;
}

sub _reference_sequence_matches {
    my $self = shift;

    my $roi_list = $self->build->region_of_interest_set;
    my $roi_reference = $roi_list->reference;
    my $reference = $self->build->reference_sequence_build;

    unless($roi_reference) {
        $self->error_message('no reference set on region of interest ' . $roi_list->name);
        return;
    }

    unless ($roi_reference->is_compatible_with($reference)) {
        if(Genome::Model::Build::ReferenceSequence::Converter->get(source_reference_build => $roi_reference, destination_reference_build => $reference)) {
            $self->status_message('Will run converter on ROI list.');
        } else {
            $self->error_message('reference sequence: ' . $reference->name . ' does not match the reference on the region of interest: ' . $roi_reference->name);
            return;
        }
    }

    return 1;
}

sub should_run {
    my $self = shift;

    return unless $self->build->region_of_interest_set;

    my $mode = $self->mode;

    my $sample_acc = $mode . '_sample';

    return $self->build->$sample_acc;

}

sub params_for_result {
    my $self = shift;
    my $build = $self->build;
    my $pp = $build->processing_profile;

    my $alignment_result;
    if($self->mode eq 'tumor') {
        $alignment_result = $build->merged_alignment_result;
    } else {
        $alignment_result = $build->control_merged_alignment_result;
    }

    unless($alignment_result) {
        die $self->error_message('No alignment result found for ' . $self->mode);
    }

    my $fl = $build->region_of_interest_set;

    my $use_short_roi = 1;
    my $short_roi_input = Genome::Model::Build::Input->get(name => 'short_roi_names', build => $build);
    if($short_roi_input) {
        $use_short_roi = $short_roi_input->value_id;
    }

    my $merge_regions = 1;
    my $merge_regions_input = Genome::Model::Build::Input->get(name => 'merge_roi_set', build => $build);
    if($merge_regions_input) {
        $merge_regions = $merge_regions_input->value_id;
    }

    return (
        alignment_result_id => $alignment_result->id,
        region_of_interest_set_id => $fl->id,
        minimum_depths => $pp->refcov_minimum_depths,
        wingspan_values => $pp->refcov_wingspan_values,
        minimum_base_quality => ($pp->refcov_minimum_base_quality || 0),
        minimum_mapping_quality => ($pp->refcov_minimum_mapping_quality || 0),
        use_short_roi_names => $pp->refcov_use_short_names,
        merge_contiguous_regions => $pp->refcov_merge_roi_regions,
        roi_track_name => ($pp->refcov_roi_track_name || undef),
        test_name => ($ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef),
    );
}

sub link_result_to_build {
    my $self = shift;
    my $result = shift;
    my $build = $self->build;

    my $link = join('/', $build->data_directory, 'coverage', $self->mode);
    my $label = join('_', 'coverage_stats', $self->mode);
    Genome::Sys->create_symlink($result->output_dir, $link);
    $result->add_user(label => $label, user => $build);

    $self->reference_coverage_result($result);

    return 1;
}

1;
