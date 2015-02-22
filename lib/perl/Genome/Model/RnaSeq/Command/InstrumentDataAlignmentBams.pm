package Genome::Model::RnaSeq::Command::InstrumentDataAlignmentBams;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::Command::InstrumentDataAlignmentBams {
    is => 'Genome::Model::Command::InstrumentData::AlignmentBams',
    has => [
        build => {
            is => 'Genome::Model::Build::RnaSeq',
            shell_args_position => 1,
        },
    ],
};

1;
