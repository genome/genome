package Genome::Model::SomaticValidation::Command::AlignReads;

use strict;
use warnings;

use Genome;

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
            default => 'apipe',
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

sub execute {
    my $self = shift;
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

    my $build_alignment_dir = join('/', $build->data_directory, 'alignments');

    for my $r (@results) {
        my @i = $r->instrument_data;
        my $sample = $i[0]->sample;
        if($sample eq $build->tumor_sample) {
            $self->merged_alignment_result_id($r->id);
            $self->merged_bam_path($r->bam_path);
            $r->add_user(label => 'merged_alignment', user => $build);
            Genome::Sys->create_symlink($r->output_dir, $build_alignment_dir . '/tumor');
        } elsif ($sample eq $build->normal_sample) {
            $self->control_merged_alignment_result_id($r->id);
            $self->control_merged_bam_path($r->bam_path);
            $r->add_user(label => 'control_merged_alignment', user => $build);
            Genome::Sys->create_symlink($r->output_dir, $build_alignment_dir . '/normal');
        } else {
            $self->warning_message('Unexpected alignment result encountered! Check samples of instrument data.');
            $r->add_user(label => 'uses', user => $build);
        }
    }

    return 1;
}

1;

