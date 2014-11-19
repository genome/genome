package Genome::Command::Base;

use strict;
use warnings;

use Genome;

class Genome::Command::Base {
    is_abstract => 1,
    is => 'Command::Tree',
};

use constant ALTERNATE_FROM_CLASS => (
    # find_class => via_class => via_class_methods
    # first method is the default method
    # the default method is used automatically if not the paramater
    # data type so it should be the most verbose option
    'Genome::Sample' => {
        'Genome::PopulationGroup' => ['samples'],
    },
    'Genome::InstrumentData' => {
        'Genome::Model' => ['instrument_data'],
        'Genome::Model::Build' => ['instrument_data'],
        'Genome::ModelGroup' => ['instrument_data'],
    },
    'Genome::Model' => {
        'Genome::Model::Build' => ['model'],
        'Genome::ModelGroup' => ['models'],
        'Genome::Config::AnalysisProject' => ['models'],
    },
    'Genome::Model::Build' => {
        'Genome::Model' => ['builds'],
    },
);

sub execute_with_shell_params_and_exit {
    my $self = shift;

    local %Command::V2::ALTERNATE_FROM_CLASS = (
        %Command::V2::ALTERNATE_FROM_CLASS,
        ALTERNATE_FROM_CLASS,
    );

    $self->SUPER::execute_with_shell_params_and_exit(@_);
}

1;
