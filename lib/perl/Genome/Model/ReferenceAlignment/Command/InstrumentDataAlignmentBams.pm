package Genome::Model::ReferenceAlignment::Command::InstrumentDataAlignmentBams;

use strict;
use warnings;

use Genome;

class Genome::Model::ReferenceAlignment::Command::InstrumentDataAlignmentBams {
    is => 'Genome::Model::Command::InstrumentData::AlignmentBams',
    has => [
        build => {
            is => 'Genome::Model::Build::ReferenceAlignment',
            shell_args_position => 1,
        },
    ],
};

1;
