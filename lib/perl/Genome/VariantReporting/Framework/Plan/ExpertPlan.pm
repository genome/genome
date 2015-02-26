package Genome::VariantReporting::Framework::Plan::ExpertPlan;

use strict;
use warnings FATAL => 'all';
use Genome;
use Genome::VariantReporting::Framework::Plan::Base;

class Genome::VariantReporting::Framework::Plan::ExpertPlan {
    is => 'Genome::VariantReporting::Framework::Plan::Base',
    has => [
        run_params => {
            is => 'HASH',
            default => {},
        },
    ],
};

sub category {
    'experts';
}

sub needed_translations {
    my $self = shift;
    my $needed = Set::Scalar->new(map {$self->run_params->{$_}} $self->get_class->run_class->required_translated_input_names);
    $needed->insert(map {@{$self->run_params->{$_}}} $self->get_class->run_class->required_translated_is_many_input_names);
    return $needed;
}

sub translate {
    my ($self, $translations) = @_;
    return $self->_translate($translations, 'run_params', $self->get_class->run_class);
}

# ExpertPlans don't have any params but have run_params instead
sub as_hashref {
    my $self = shift;

    my %body;
    while (my ($name, $value) = each %{$self->run_params}) {
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

    eval {
        $self->get_class->run_class;
    };

    if (my $error = $@) {
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => [],
            desc => $error,
        );
    }

    return @errors;
}

sub __object_errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__object_errors__;
    push @errors, $self->get_class->run_class->__planned_errors__($self->run_params);
    return @errors;
}

1;
