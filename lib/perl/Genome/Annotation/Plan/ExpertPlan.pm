package Genome::Annotation::Plan::ExpertPlan;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Plan::ExpertPlan {
    is => 'Genome::Annotation::Plan::Base',
};

sub category {
    'experts';
}

sub object {
    my $self = shift;
    return $self->factory->get_object($self->category,
        $self->name, {});
}

sub adaptor_object {
    my $self = shift;
    my $adaptor_class = $self->object->adaptor_class;
    return $adaptor_class->create();
}

sub validate_object {
    my $self = shift;
    $self->adaptor_object->validate_with_plan_params($self->params);
}

1;
