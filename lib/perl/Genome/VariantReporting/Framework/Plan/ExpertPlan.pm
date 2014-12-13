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

sub needed_translations {
    my $self = shift;
    return Set::Scalar->new(map {$self->adaptor_params->{$_}}
        $self->get_class->adaptor_class->required_translated_input_names);
}

sub translate {
    my ($self, $translations) = @_;
    return $self->_translate($translations, 'adaptor_params', $self->get_class->adaptor_class);
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
    push @errors, $self->get_class->adaptor_class->__planned_output_errors__($self->adaptor_params);
    return @errors;
}

1;
