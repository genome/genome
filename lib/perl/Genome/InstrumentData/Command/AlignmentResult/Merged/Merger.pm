package Genome::InstrumentData::Command::AlignmentResult::Merged::Merger;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::AlignmentResult::Merged::Merger {
    is => 'Command',

    has_input => [
        input_bams => {
            is => 'Text',
            is_many => 1,
            doc => 'paths to the BAM files to be merged',
        },
        output_path => {
            is => 'Text',
            is_output => 1,
            doc => 'path to write the merged BAM',
        },
        parameters => {
            is => 'Text',
            doc => 'string representing parameters for the merger',
            is_optional => 1,
        },
        version => {
            is => 'Text',
            doc => 'version of the merger to use',
        },
        scratch_directory => {
            is => 'Text',
            doc => 'path to a directory for scratch work to be done',
            is_optional => 1,
        },
        samtools_version => {
            is => 'Text',
            doc => 'version of samtools to use',
        }
    ],
};

1;
