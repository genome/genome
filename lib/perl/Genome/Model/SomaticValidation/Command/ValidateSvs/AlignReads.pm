package Genome::Model::SomaticValidation::Command::ValidateSvs::AlignReads;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticValidation::Command::ValidateSvs::AlignReads {
    is => 'Command::V2',
    has_input => [
        build_id => {
            is => 'Number',
            doc => 'id of build for which to run alignments',
            is_output => 1,
        },
        reference_build_id => {
            is => 'Text',
            doc => 'id of the reference build against which to run alignments',
            is_optional => 1,
        },
        skip => {
            is => 'Boolean',
            doc => 'signal from first step whether or not to run',
            default_value => 0,
        },
    ],
    has => [
        build => {
            is => 'Genome::Model::Build',
            id_by => 'build_id',
            doc => 'build for which to run alignments',
        },
        reference_build => {
            is => 'Genome::Model::Build',
            id_by => 'reference_build_id',
            is_optional => 1,
        },
        output_dir => {
            is_output => 1,
            is_input => 1,
            is => 'Text',
            doc => 'Place where the output goes',
        },

    ],
    has_param => [
        lsf_queue => {
            default => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT},
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

    if($self->skip) {
        $self->status_message("skip signal received. not running.");
        return 1;
    } elsif(not $self->reference_build) {
        die $self->error_message("A reference is required when not skipping this step, but one was not found.");
    }



    my @instrument_data = $build->instrument_data;
    my $result = Genome::InstrumentData::Composite->get_or_create(
        inputs => {
            instrument_data => \@instrument_data,
            reference_sequence_build => $self->reference_build,
        },
        strategy => 'instrument_data aligned to reference_sequence_build using bwa 0.5.9 [-t 4 -q 5::] then merged using picard 1.46 then deduplicated using picard 1.46 api v2',#TODO make me a processing profile parameter, backfilled default
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

    my $build_alignment_dir = join('/', $self->output_dir, 'alignments');
    Genome::Sys->create_directory($build_alignment_dir);

    for my $r (@results) {
        my @i = $r->instrument_data;
        my $sample = $i[0]->sample;
        if($sample eq $build->tumor_sample) {
            $self->merged_alignment_result_id($r->id);
            $self->merged_bam_path($r->merged_alignment_bam_path);
            $r->add_user(label => 'sv_validation_merged_alignment', user => $build);
            Genome::Sys->create_symlink($r->output_dir, $build_alignment_dir . '/tumor');
        } elsif ($sample eq $build->normal_sample) {
            $self->control_merged_alignment_result_id($r->id);
            $self->control_merged_bam_path($r->merged_alignment_bam_path);
            $r->add_user(label => 'sv_validation_control_merged_alignment', user => $build);
            Genome::Sys->create_symlink($r->output_dir, $build_alignment_dir . '/normal');
        } else {
            $self->warning_message('Unexpected alignment result encountered! Check samples of instrument data.');
            $r->add_user(label => 'uses', user => $build);
        }
    }

    return 1;
}

1;

