package Genome::Model::Build::SomaticInterface;

use strict;
use warnings;
use Genome;

class Genome::Model::Build::SomaticInterface {
    is => 'UR::Object',
    is_abstract => 1,
};

sub reference_sequence_build {
    my $self = shift;
    $self->fatal_message('Abstract: (reference_sequence_build) needs to be defined on class (%s)', $self->class);
}

sub individual_common_name {
    my $self = shift;
    $self->fatal_message('Abstract: (individual_common_name) needs to be defined on class (%s)', $self->class);
}

1;
