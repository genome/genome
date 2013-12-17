package Genome::Model::RnaSeq::Command::AlignReads;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::Command::AlignReads {
    is => 'Command::V2',
    has_input => [ build_id => { is => 'Number', doc => 'id of build for which to run alignments',
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
            default => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT},
        },
    ],
    has_transient_optional_output => [
        merged_alignment_result_id => {
            is => 'Number',
            doc => 'id of the merged alignment result for the instrument data',
        },
        merged_bam_path => {
            is => 'Text',
            doc => 'Path to the merged instrument data bam',
        },
        individual_alignment_results => {
            is => 'Genome::InstrumentData::AlignmentResult',
            is_many => 1,
            doc => 'original alignment results per-instrument-data',
        },
    ],
    doc => 'align reads',
};

sub sub_command_category { 'pipeline' }

sub execute {
    my $self = shift;
    my $build = $self->build;

    my @instrument_data = $build->instrument_data;
    my $result = Genome::InstrumentData::Composite->get_or_create(
        inputs => {
            instrument_data => \@instrument_data,
            reference_sequence_build => $build->reference_sequence_build,
            picard_version => $build->processing_profile->picard_version,
            trimmer_name => $build->processing_profile->read_trimmer_name,
            trimmer_version => $build->processing_profile->read_trimmer_version,
            trimmer_params => $build->processing_profile->read_trimmer_params,
            annotation_build => $build->annotation_build,
        },
        strategy => $self->_generate_alignment_strategy,
        log_directory => $build->log_directory,
    );

    my @bams = $result->bam_paths;

    unless(scalar(@bams) == 1) {
        die($self->error_message('Found ' . scalar(@bams) . ' from alignment when 1 was expected. This model will probably fail.'));
    }

    $self->status_message("Alignment BAM paths:\n " . join("\n ", @bams));

    my @results = $result->_merged_results;

    my $build_alignment_dir = join('/', $build->data_directory, 'alignments');

    for my $r (@results) {
        $self->merged_alignment_result_id($r->id);
        $self->merged_bam_path($r->merged_alignment_bam_path);
        $r->add_user(label => 'merged_alignment', user => $build);
        Genome::Sys->create_symlink($r->output_dir, $build_alignment_dir);

        my @individual_alignment_results = $r->collect_individual_alignments;
        for my $i (@individual_alignment_results) {
            $r->add_user(label => 'individual_alignment', user => $build);
        }
        $self->individual_alignment_results(\@individual_alignment_results);
    }

    return 1;
}

sub _generate_alignment_strategy {
    my $self = shift;
    my $pp = $self->build->processing_profile;
    my $picard_version = $pp->picard_version;
    my $aligner = 'per-lane-tophat';
    my $aligner_version = $pp->read_aligner_version;
    my $aligner_params = $pp->read_aligner_params;
    my $strategy = "instrument_data aligned to reference_sequence_build and annotation_build using $aligner $aligner_version [$aligner_params] then merged using picard $picard_version";
    if (defined($pp->deduplication_handler)) {
        if ($pp->deduplication_handler eq 'picard') {
            $strategy .= ' then deduplicated using '. $pp->deduplication_handler .' '. $picard_version;
        } else {
            die('Failed to recognize the deduplication_handler '. $pp->deduplication_handler);
        }
    }
    return $strategy .' api v2';
}

sub max_elapsed_log_time { 7 * 24 * 3600 };

1;

