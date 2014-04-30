package Genome::Annotation::Filter::Base;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Filter::Base {
    is => 'Genome::Annotation::ComponentBase',
    is_abstract => 1,
    attributes_have => [
        is_translated => {is => 'Boolean', default => 0},
    ],
};

sub requires_experts {
    return ();
}

sub translated_inputs {
    my $self = shift;
    my $meta = $self->__meta__;
    return map {$_->property_name} $meta->properties(
        is_translated => 1,
    );
}

1;
