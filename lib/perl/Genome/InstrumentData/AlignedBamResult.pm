package Genome::InstrumentData::AlignedBamResult;

use Genome;

use warnings;
use strict;

class Genome::InstrumentData::AlignedBamResult {
    is_abstract => 1,
    is => 'Genome::SoftwareResult',

    has_output => [
        bam_file => {
            is => "Text",
            doc => "The path to the aligned bam file.",
        },
    ]
};


1;

