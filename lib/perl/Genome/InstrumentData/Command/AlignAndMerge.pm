package Genome::InstrumentData::Command::AlignAndMerge;

use strict;
use warnings;
use Genome;

class Genome::InstrumentData::Command::AlignAndMerge {
    is => ['Command::V2'],
    has => [
        name => {
            is => 'Text',
            doc => 'The aligner to use',
        },
        version => {
            is => 'Text',
            doc => 'Version of the aligner to use',
        },
        params => {
            is => 'Text',
            doc => 'Parameters to pass to the aligner',
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
            doc => 'Id of the data to align',
        },
        reference_sequence_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'Id of the reference to which to align',
        },
        picard_version => {
            is => 'Text',
            doc => 'The version of Picard to use when needed by aligners/filters',
        },
        samtools_version => {
            is => 'Text',
            doc => 'The version of Samtools to use when needed by aligners/filters',
        },
        result_users => {
            is => 'HASH',
            doc => 'mapping of labels to user objects. Will be added to any generated results',
        },
    ],
    has_optional_input => [
        bedtools_version => {
            is => 'Text',
            doc => 'The version of Bedtools to use when needed by aligners/filters',
        },
        annotation_build_id => {
            is => 'Number',
            doc => 'Id of the annotation build to use when aligning',
            is_optional => 1,
        },
        instrument_data_segment_id => {
            is => 'Text',
            doc => 'A specific segment of the data to align',
        },
        instrument_data_segment_type => {
            is => 'Text',
            doc => 'How the data is segmented',
        },
        instrument_data_filter => {
            is => 'Text',
            valid_values => [ 'forward-only', 'reverse-only', undef ],
        },
        force_fragment => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Treat all reads as fragment reads',
        },
        trimmer_name => {
          is => 'Text',
          doc => 'name of the read trimmer to use before alignment is performed'
        },
        trimmer_params => {
          is => 'Text',
          doc => 'params to use with the read trimmmer'
        },
        trimmer_version => {
          is => 'Text',
          doc => 'version of the read trimmer to use'
        },
    ],
    has_optional_output => [
        result_id => {
            is => 'Text',
            doc => 'The ID of the result generated/found when running the command',
        },
        alignment_result => {
            is => 'Genome::InstrumentData::AlignmentResult::Merged',
            id_by => 'result_id',
            doc => 'The result generated/found when running the command',
        }
    ],
};

sub execute {
    my $self = shift;

    my $result = $self->_process_alignments('get_or_create');
    unless ($result) {
        $self->error_message("Error finding or generating alignments!");
        return 0;
    }
    $self->result_id($result->id);

    my $per_lane_results = $self->_process_per_lane_alignments('get_or_create');
    unless ($per_lane_results) {
        $self->error_message("Error finding or generating per-lane alignments!");
        return 0;
    }

    return $result;
}

sub shortcut {
    my $self = shift;

    my $result = $self->_process_alignments('get_with_lock');
    unless ($result) {
        return undef;
    }
    $self->result_id($result->id);

    my $per_lane_results = $self->_process_per_lane_alignments('get_with_lock');
    unless ($per_lane_results) {
        return undef;
    }

    return $result;
}

sub _process_alignments {
    my $self = shift;
    my $mode = shift;

    my %params = $self->_alignment_params;
    my $result = Genome::InstrumentData::AlignmentResult::Merged::Speedseq->$mode(
        instrument_data => [$self->instrument_data],
        %params,
    );

    return $result;
}

sub _process_per_lane_alignments {
    my $self = shift;
    my $mode = shift;

    my %params = $self->_alignment_params;
    my @per_lane_results;
    for my $instrument_data ($self->instrument_data) {
        push @per_lane_results, Genome::InstrumentData::AlignmentResult::Speedseq->$mode(
            instrument_data => $instrument_data,
            samtools_version => $self->samtools_version,
            picard_version => $self->picard_version,
            %params,
        );
    }
    return (scalar(@per_lane_results) == scalar($self->instrument_data));
}

sub _alignment_params {
    my $self = shift;

    return (
        reference_build => $self->reference_sequence_build,
        users => $self->result_users,
        aligner_name => $self->name,
        aligner_version => $self->version,
        aligner_params => $self->params,
    );
}

1;
