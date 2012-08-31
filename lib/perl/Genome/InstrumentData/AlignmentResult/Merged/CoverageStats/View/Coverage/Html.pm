package Genome::InstrumentData::AlignmentResult::Merged::CoverageStats::View::Coverage::Html;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::AlignmentResult::Merged::CoverageStats::View::Coverage::Html {
    is => 'Genome::View::Status::Html',
    has => [
        perspective => { value => 'coverage' },
    ],
};

1;
