package Genome::Model::MutationalSignificance::Command::CreateClinicalData;

use strict;
use warnings;

use Genome;

class Genome::Model::MutationalSignificance::Command::CreateClinicalData {
    is => ['Command::V2'],
    has_input_output => [
        clinical_data_file => {
            is => 'String',
        },
    ],
};

sub execute {
    my $self = shift;

    $self->clinical_data_file("a_clinical_data_file");
    $self->status_message('Created clinical data file');
    return 1;
}

1;
