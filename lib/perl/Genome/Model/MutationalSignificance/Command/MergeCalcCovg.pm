package Genome::Model::MutationalSignificance::Command::MergeCalcCovg;

use strict;
use warnings;

use Genome;

class Genome::Model::MutationalSignificance::Command::MergeCalcCovg {
    is => ['Command::V2'],
    has_input => [
        output_files => {
            is => 'Text',
            is_many => 1,
        },
    ],
    has_input_output => [
        output_dir => {
            is => 'Text',
        },
    ],
};

#Doesn't do anything except confirm that all of the CalcCovg operations completed
sub execute {
    my $self = shift;

    my $output_file_count = scalar $self->output_files;

    $self->status_message("MergeCalcCovg detected $output_file_count output files");

    return 1;
}

1;
