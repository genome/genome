package Genome::VariantReporting::Framework::Component::Reporter;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Framework::Component::Reporter {
    is => 'Genome::VariantReporting::Framework::Component::Base',
    is_abstract => 1,
};

sub name {
    die "abstract";
}

sub requires_interpreters {
    die "abstract - must return a list of one or more interpreter names";
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
