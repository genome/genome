package Genome::Model::SomaticValidation::Command::AlignReads;

use strict;
use warnings;

use Genome;

use File::Spec;

class Genome::Model::SomaticValidation::Command::AlignReads {
    is => 'Command::V2',
    has_input => [
        build_id => {
            is => 'Number',
            doc => 'id of build for which to run alignments',
            is_output => 1,
        },
    ],
    has => [
        build => {
            is => 'Genome::Model::Build',
            id_by => 'build_id',
            doc => 'build for which to run alignments',
        },
    ],
    has_param => [
        lsf_queue => {
            default => Genome::Config::get('lsf_queue_build_worker'),
        },
    ],
    has_transient_optional_output => [
        merged_alignment_result_id => {
            is => 'Number',
            doc => 'id of the merged alignment result for the instrument data',
        },
        control_merged_alignment_result_id => {
            is => 'Number',
            doc => 'id of the merged alignment result for the control instrument data',
        },
        merged_bam_path => {
            is => 'Text',
            doc => 'Path to the merged instrument data bam',
        },
        control_merged_bam_path => {
            is => 'Text',
            doc => 'Path to the merged control instrument data bam',
        },
    ],
    doc => 'align reads',
};

sub sub_command_category { 'pipeline steps' }

sub shortcut {
    my $self = shift;

    return $self->_find_existing_alignments;
}


sub execute {
    my $self = shift;

    return 1 if $self->_find_existing_alignments;

    my $build = $self->build;

    my @known_sites_inputs = $build->inputs(name => 'known_sites');

    my @instrument_data = $build->instrument_data;
    my $result = Genome::InstrumentData::Composite->get_or_create(
        inputs => {
            instrument_data => \@instrument_data,
            reference_sequence_build => $build->reference_sequence_build,
            (@known_sites_inputs ? (known_sites => [ map { $_->value } @known_sites_inputs ]) : ()),
        },
        strategy => $build->processing_profile->alignment_strategy,
        log_directory => $build->log_directory,
        result_users => Genome::SoftwareResult::User->user_hash_for_build($build),
    );

    my @bams = $result->bam_paths;

    my $num_expected_samples = 0;
    $num_expected_samples++ if $build->tumor_sample;
    $num_expected_samples++ if $build->normal_sample;

    unless(scalar(@bams) == $num_expected_samples) {
        $self->warning_message('Found ' . scalar(@bams) . ' from alignment when ' . $num_expected_samples . ' expected. This model will probably fail.');
    }

    $self->status_message("Alignment BAM paths:\n " . join("\n ", @bams));

    my @results = $result->_merged_results;

    for my $r (@results) {
        my @i = $r->instrument_data;
        my $sample = $i[0]->sample;
        if($sample eq $build->tumor_sample) {
            $self->_assign_tumor_sample_alignment($r);
        } elsif ($sample eq $build->normal_sample) {
            $self->_assign_normal_sample_alignment($r);
        } else {
            $self->warning_message('Unexpected alignment result encountered! Check samples of instrument data.');
            $r->add_user(label => 'uses', user => $build);
        }
    }

    return 1;
}

sub _find_existing_alignments {
    my $self = shift;
    my $build = $self->build;

    my @existing_alignments = $build->prealigned_data;

    return unless @existing_alignments;

    for my $e (@existing_alignments) {
        unless ($e->isa('Genome::InstrumentData::AlignedBamResult::Merged')) {
            $self->fatal_message('Unexpected result type in prealigned data: %s', $e->class);
        }

        my $sample = $e->sample_name;
        if ($sample eq $build->tumor_sample->name) {
            $self->_assign_tumor_sample_alignment($e);
        } elsif ($sample eq $build->normal_sample->name) {
            $self->_assign_normal_sample_alignment($e);
        } else {
            $self->fatal_message('Found prealigned data that does not match build samples!: %s', $e->id);
        }

        $e->add_user(label => 'shortcut', user => $build);
        if (my $anp = $build->model->analysis_project) {
            $e->add_user(label => 'sponsor', user => $anp);
        }
    }

    return @existing_alignments;
}

sub _assign_tumor_sample_alignment {
    my $self = shift;
    my $alignment = shift;

    my $build = $self->build;
    my $build_alignment_dir = File::Spec->join($build->data_directory, 'alignments');

    $self->merged_alignment_result_id($alignment->id);
    $self->merged_bam_path($alignment->bam_path);
    $alignment->add_user(label => 'merged_alignment', user => $build);
    Genome::Sys->create_symlink($alignment->output_dir, File::Spec->join($build_alignment_dir, 'tumor'));

    return 1;
}

sub _assign_normal_sample_alignment {
    my $self = shift;
    my $alignment = shift;

    my $build = $self->build;
    my $build_alignment_dir = File::Spec->join($build->data_directory, 'alignments');

    $self->control_merged_alignment_result_id($alignment->id);
    $self->control_merged_bam_path($alignment->bam_path);
    $alignment->add_user(label => 'control_merged_alignment', user => $build);
    Genome::Sys->create_symlink($alignment->output_dir, File::Spec->join($build_alignment_dir, 'normal'));

    return 1;
}

1;

