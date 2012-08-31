package Genome::InstrumentData::IntermediateAlignmentResult::Command::Bwa;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::IntermediateAlignmentResult::Command::Bwa {
    is => 'Command::V2',
    has_input => [
        instrument_data_id => {
            is => 'Text',
            doc => 'The instrument data upon which to run `bwa aln`',
        },
        aligner_name => {
            is => 'Text',
            doc => 'The aligner to use (must be "bwa")',
            valid_values => ['bwa'],
            default_value => 'bwa',
        },
        aligner_version => {
            is => 'Text',
            doc => 'Version of bwa to use',
        },
        aligner_params => {
            is => 'Text',
            doc => 'Additional parameters to pass to `bwa aln`',
        },
        aligner_index_id => {
            is => 'Text',
            doc => 'The index for the reference',
        },
        input_file => {
            is => 'Text',
            doc => 'The file to align',
        },
        input_pass => {
            is => 'Number',
            doc => 'Direction of reads to align',
        },
        instrument_data_segment_type => {
            is => 'Text',
            doc => 'Parameter by which to split a large set of reads',
        },
        instrument_data_segment_id => {
            is => 'Text',
            doc => 'The specific value of the segment type to work with in this instance',
        },
        samtools_version => {
            is => 'Text',
            doc => 'The samtools version used',
        },
        test_name => {
            is => 'Text',
            doc => 'A test name for the software result',
            default_value => $ENV{SOFTWARE_RESULT_TEST_NAME},
        },
    ],
    has_optional_output => [
        _result => {
            is => 'Genome::InstrumentData::IntermediateAlignmentResult::Bwa',
            doc => 'the software result gotten or created',
        }
    ],
};

sub execute {
    my $self = shift;

    my %intermediate_params = (
        instrument_data_id => $self->instrument_data_id,
        aligner_name => $self->aligner_name,
        aligner_version => $self->aligner_version,
        aligner_params => $self->aligner_params,
        aligner_index_id => $self->aligner_index_id,
        input_file => $self->input_file,
        input_pass => $self->input_pass,
        instrument_data_segment_type => $self->instrument_data_segment_type,
        instrument_data_segment_id => $self->instrument_data_segment_id,
        samtools_version => $self->samtools_version,
        test_name => $self->test_name,
    );

    my $result = Genome::InstrumentData::IntermediateAlignmentResult::Bwa->get_or_create(
        %intermediate_params
    );

    unless($result) {
        die $self->error_message("Failed to generate IntermediateAlignmentResult. Params were: " . Data::Dumper::Dumper(\%intermediate_params));
    }

    $self->_result($result);

    return 1;
}

1;
