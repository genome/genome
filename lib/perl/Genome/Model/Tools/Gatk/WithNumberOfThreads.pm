package Genome::Model::Tools::Gatk::WithNumberOfThreads;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Gatk::WithNumberOfThreads {
    is => 'UR::Object',
    is_abstract => 1,
    attributes_have => {
        is_input => { is => 'Boolean', is_optional => 1 },        
    },
    has_optional_input => {
        number_of_threads => {
            is => 'Number',
            doc => 'Controls the number of data threads sent to the processor',
        },
    },
};

sub number_of_threads_param_for_java_command {
    my $self = shift;

    return '' if not $self->number_of_threads;

    return ' -nt '.$self->number_of_threads;
}

1;

