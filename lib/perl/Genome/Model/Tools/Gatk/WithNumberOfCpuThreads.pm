package Genome::Model::Tools::Gatk::WithNumberOfCpuThreads;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Gatk::WithNumberOfCpuThreads {
    is => 'UR::Object',
    is_abstract => 1,
    attributes_have => {
        is_input => { is => 'Boolean', is_optional => 1 },        
    },
    has_optional_input => {
        number_of_cpu_threads => {
            is => 'Number',
            doc => 'Controls the number of CPU threads allocated to each data thread',
        },
    },
};

sub number_of_cpu_threads_param_for_java_command {
    my $self = shift;

    return '' if not $self->number_of_cpu_threads;

    return ' -nct '.$self->number_of_cpu_threads;
}

1;

