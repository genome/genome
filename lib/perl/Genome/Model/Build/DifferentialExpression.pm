package Genome::Model::Build::DifferentialExpression;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::DifferentialExpression {
    is => 'Genome::Model::Build',
};

sub transcript_convergence_directory {
    my $self = shift;
    return $self->data_directory . '/transcript_convergence';
}

1;

