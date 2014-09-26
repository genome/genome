package Genome::VariantReporting::Command::CombineSnvIndelReports;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Command::CombineSnvIndelReports {
    is => 'Command::V2',
    has_input => [
        input_directories => {
            is_many => 1,
            is => 'Path',
        },
        file_name => {
            is => 'Text',
        },
    ],
    has_calculated => [
        reports => {
            is_output => 1,
            is_many => 1,
            is => 'File',
            calculate => q|map {File::Spec->join($_, $self->file_name)} $self->input_directories|,
        },
    ],
};

1;

