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
            is => 'Genome::InstrumentData::AlignedBamResult::Merged',
            id_by => 'result_id',
            doc => 'The result generated/found when running the command',
        },
        per_lane_alignment_result_ids => {
            is => 'Text',
            doc => 'The IDs of the per-lane results generated/found when running the command',
            is_many => 1,
        },
        per_lane_alignment_results => {
            is => 'Genome::InstrumentData::AlignmentResult',
            is_many => 1,
            calculate => q{
                return Genome::InstrumentData::AlignmentResult->get([$self->per_lane_alignment_result_ids]);
            },
            doc => 'The per-lane results generated/found when running the command',
        },
    ],
};

sub execute {
    my $self = shift;

    for my $param (qw(trimmer_name trimmer_params trimmer_version)) {
        if (defined($self->$param)) {
            $self->error_message('Optional parameter (%s) not supported at this time.', $param);
            return 0;
        }
    }

    my $result = $self->_process_alignments('get_or_create');
    unless ($result) {
        $self->error_message("Error finding or generating alignments!");
        return 0;
    }
    $self->result_id($result->id);

    my @per_lane_results = $self->_process_per_lane_alignments('get_or_create');
    unless (scalar(@per_lane_results) == scalar($self->instrument_data)) {
        $self->error_message("Error finding or generating per-lane alignments!");
        return 0;
    }
    $self->per_lane_alignment_result_ids([map { $_->id } @per_lane_results]);

    return $result;
}

sub shortcut {
    my $self = shift;

    my $result = $self->_process_alignments('get_with_lock');
    unless ($result) {
        return undef;
    }
    $self->result_id($result->id);

    my @per_lane_results = $self->_process_per_lane_alignments('get_with_lock');
    unless (scalar(@per_lane_results) == scalar($self->instrument_data)) {
        return undef;
    }

    return $result;
}

sub _process_alignments {
    my $self = shift;
    my $mode = shift;

    my $class = $self->merged_result_class;
    my %params = $self->_alignment_params;
    my $result = $class->$mode(
        instrument_data => [$self->instrument_data],
        %params,
    );

    return $result;
}

sub _process_per_lane_alignments {
    my $self = shift;
    my $mode = shift;

    my $class = $self->per_lane_result_class;
    my %params = $self->_alignment_params;
    my @per_lane_results;
    for my $instrument_data ($self->instrument_data) {
        push @per_lane_results, $class->$mode(
            instrument_data => $instrument_data,
            samtools_version => $self->samtools_version,
            picard_version => $self->picard_version,
            %params,
        );
    }

    return @per_lane_results;
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

sub merged_result_class {
    my $self = shift;
    return 'Genome::InstrumentData::AlignmentResult::Merged::' . ucfirst($self->name);
}

sub per_lane_result_class {
    my $self = shift;
    return 'Genome::InstrumentData::AlignmentResult::' . ucfirst($self->name);
}

1;
