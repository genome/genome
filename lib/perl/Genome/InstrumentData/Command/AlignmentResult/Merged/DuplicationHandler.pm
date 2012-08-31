package Genome::InstrumentData::Command::AlignmentResult::Merged::DuplicationHandler;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::AlignmentResult::Merged::DuplicationHandler {
    is => 'Command',

    has_input => [
        input_bam => {
            is => 'Text',
            doc => 'path to the BAM file in which duplicates should be handled',
        },
        output_path => {
            is => 'Text',
            is_output => 1,
            doc => 'path to write the BAM with duplications handled',
        },
        parameters => {
            is => 'Text',
            doc => 'string representing parameters for the duplication handler',
            is_optional => 1,
        },
        version => {
            is => 'Text',
            doc => 'version of the duplication handler to use',
        },
        scratch_directory => {
            is => 'Text',
            doc => 'path to a directory for scratch work to be done',
            is_optional => 1,
        },
        log_file => {
            is => 'Text',
            doc => 'path to direct output from duplication handler',
        },
        metrics_file => {
            is => 'Text',
            doc => 'path to location to output metrics file',
        },
    ],
};

1;
