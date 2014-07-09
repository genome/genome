package Genome::VariantReporting::Framework::Component::Reporter;

use strict;
use warnings FATAL => 'all';
use Genome;
use Carp qw(confess);

class Genome::VariantReporting::Framework::Component::Reporter {
    is => 'Genome::VariantReporting::Framework::Component::Base',
    is_abstract => 1,
};

sub name {
    my $self = shift;
    my $class = $self->class;
    confess "Abstract method 'name' must be defined in class '$class'";
}

sub requires_interpreters {
    my $self = shift;
    my $class = $self->class;
    confess "Abstract method 'requires_interpreters' must be defined in class '$class'";
}

sub initialize {
    # this gets called before interpretations are given to ->report method
    return;
}

sub report {
    my ($self, $interpretations) = shift;
    # do something with the interpretations to produce one or more reports
    return;
}

sub finalize {
    # this gets called after interpretations are given to ->report method
    return;
}

1;
