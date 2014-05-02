package Genome::Annotation::Interpreter::Base;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Interpreter::Base {
    is => 'Genome::Annotation::ComponentBase',
    is_abstract => 1,
    attributes_have => [
        is_translated => {is => 'Boolean', default => 0},
    ],
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

sub translated_inputs {
    my $self = shift;
    my $meta = $self->__meta__;
    return map {$_->property_name} $meta->properties(
        is_translated => 1,
    );
}


1;
