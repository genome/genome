package Genome::Model::ReferenceAlignment::Command::AlignReads;

use strict;
use warnings;

use File::Spec;

use Genome;

class Genome::Model::ReferenceAlignment::Command::AlignReads {
    is => 'Genome::Model::ReferenceAlignment::Command::PipelineBase',
    has_param => [
        lsf_queue => {
            default => Genome::Config::get('lsf_queue_build_worker'),
        },
    ],
    has_transient_optional_output => [
        merged_alignment_result => {
            is => 'Genome::InstrumentData::AlignedBamResult::Merged',
            doc => 'id of the merged alignment result for the instrument data',
        },
        merged_bam_path => {
            is => 'Text',
            doc => 'Path to the merged instrument data bam',
        },
        individual_alignment_results => {
            is => 'Genome::InstrumentData::AlignedBamResult',
            doc => 'the individual (usually "per-lane") alignment results for the instrument data',
            is_many => 1,
        },

    ],
    doc => 'runs the alignment dispatcher on the instrument data for the build',
};

sub execute {
    my $self = shift;
    my $build = $self->build;

    my $inputs = $self->_inputs_for_alignment;
    my $strategy = $self->_generate_strategy;
    my $composite = Genome::InstrumentData::Composite->get_or_create(
        inputs => $inputs,
        strategy => $strategy,
        log_directory => $build->log_directory,
        result_users => Genome::SoftwareResult::User->user_hash_for_build($build),
        merge_group => 'all', #RefAlign should always produce one BAM regardless of data assignment
    );

    my @bams = $composite->bam_paths;

    unless(scalar(@bams) == 1) {
        $self->fatal_message('Found %s from alignment when one expected.', scalar(@bams));
    }

    $self->merged_bam_path($bams[0]);
    $self->status_message("Alignment BAM path:\n " . $bams[0]);

    my ($result) = $composite->_merged_results;
    $self->merged_alignment_result($result);

    my $build_alignment_dir = File::Spec->join($build->data_directory, 'alignments');
    Genome::Sys->create_symlink($result->output_dir, $build_alignment_dir);
    $result->add_user(label => 'merged_alignment', user => $build);

    my @individual_alignments = $result->collect_individual_alignments;
    $self->individual_alignment_results(\@individual_alignments);

    return 1;
}

sub _inputs_for_alignment {
    my $self = shift;
    my $build = $self->build;

    my @instrument_data = $build->instrument_data;
    my $pp = $build->processing_profile;

    return {
        instrument_data => \@instrument_data,
        reference_sequence_build => $build->reference_sequence_build,
        force_fragment => $pp->force_fragment || undef,
        trimmer_name => $pp->read_trimmer_name || undef,
        trimmer_version => $pp->read_trimmer_version || undef,
        trimmer_params => $pp->read_trimmer_params || undef,
        picard_version => $pp->picard_version || undef,
        samtools_version => $pp->samtools_version || undef,
        bedtools_version => $pp->bedtools_version || undef,
    };
}

sub _generate_strategy {
    my $self = shift;
    my $build = $self->build;

    my $pp = $build->processing_profile;

    my $strategy = sprintf(
        'instrument_data aligned to reference_sequence_build using %s %s [%s] ',
        $pp->read_aligner_name,
        $pp->read_aligner_version,
        $pp->read_aligner_params,
    );

    if ($pp->merger_name) {
        $strategy .= sprintf(
            'then merged using %s %s [%s] ',
            $pp->merger_name,
            $pp->merger_version,
            $pp->merger_params,
        );
    }

    if ($pp->duplication_handler_name) {
        $strategy .= sprintf(
            'then deduplicated using %s %s [%s] ',
            $pp->duplication_handler_name,
            $pp->duplication_handler_version,
            $pp->duplication_handler_params,
        );
    }

    $strategy .= 'api v2'; #required by strategy, but overridden by values set in _inputs_for_alignment

    return $strategy;
}

1;
