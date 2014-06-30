package Genome::VariantReporting::Framework::Component::Interpreter;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Framework::Component::Interpreter {
    is => ['Genome::VariantReporting::Framework::Component::Base', 'Genome::VariantReporting::Framework::Component::WithTranslatedInputs'],
    is_abstract => 1,
};

sub name {
    my $self = shift;
    my $class = $self->class;
    die "Abstract method 'name' must be defined in class '$class'";
}

sub requires_experts {
    my $self = shift;
    my $class = $self->class;
    die "Abstract method 'requires_experts' must be defined in class '$class'";
}

sub interpret_entry {
    my $self = shift;
    my $class = $self->class;
    die "Abstract method 'interpret_entry' must be defined in class '$class'";
}

1;
