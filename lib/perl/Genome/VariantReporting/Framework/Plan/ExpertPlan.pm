package Genome::VariantReporting::Framework::Plan::ExpertPlan;

use strict;
use warnings FATAL => 'all';
use Genome;
use Genome::VariantReporting::Framework::Plan::Base;

class Genome::VariantReporting::Framework::Plan::ExpertPlan {
    is => 'Genome::VariantReporting::Framework::Plan::Base',
    has => [
        adaptor_params => {
            is => 'HASH',
            default => {},
        },
    ],
};

sub category {
    'experts';
}

sub resources_required {
    my $self = shift;

    return $self->get_class->resources_required;
}

sub adaptor_object {
    my $self = shift;
    my $adaptor_class = $self->object->adaptor_class;
    return $adaptor_class->create();
}

# ExpertPlans don't have any params but have adaptor_params instead
sub as_hashref {
    my $self = shift;

    my %body;
    while (my ($name, $value) = each %{$self->adaptor_params}) {
        $body{$name} = $value;
    }

    my %result;
    $result{$self->name} = \%body;

    return \%result;
}

sub __class_errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__class_errors__;
    return @errors if @errors; # Can't know anything else in this case

    for my $accessor ('run_class', 'adaptor_class') {
        eval {
            $self->get_class->$accessor;
        };

        if (my $error = $@) {
            push @errors, UR::Object::Tag->create(
                type => 'error',
                properties => [],
                desc => $error,
            );
        }
    }

    return @errors;
}

sub __object_errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__object_errors__;
    push @errors, $self->object->adaptor_class->__planned_output_errors__($self->adaptor_params);
    return @errors;
}

sub __provider_errors__ {
    my $self = shift;
    my $provider = shift;

    return $self->object->adaptor_class->__provided_output_errors__(
        $provider->attributes);
}

1;
