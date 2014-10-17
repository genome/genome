package Genome::VariantReporting::Framework::Component::Adaptor;

use strict;
use warnings;
use Genome;
use Params::Validate qw(validate validate_pos :types);
use Exception::Class;

class Genome::VariantReporting::Framework::Component::Adaptor {
    is => ['Command::V2', 'Genome::VariantReporting::Framework::Component::Base', 'Genome::VariantReporting::Framework::Component::WithTranslatedInputs'],
    is_abstract => 1,
    attributes_have => {
        is_planned => {
            is => "Boolean",
            default => 0,
        },
        is_provided => {
            is => "Boolean",
            default => 0,
        },
    },
    has_input => [
        provider_json => {
            is => 'Text',
        },
        variant_type => {
            is => 'Text',
            is_output => 1,
            valid_values => ['snvs', 'indels'],
        },
        plan_json => {
            is => 'Text',
        },
    ],
};

sub name {
    die "Abstract";
}

sub shortcut {
    my $self = shift;
    return $self->execute();
}

sub execute {
    my $self = shift;
    $self->debug_message("Resolving plan attributes");
    $self->resolve_plan_attributes;

    $self->debug_message("Resolving attributes from the resource-provider");
    $self->resolve_provided_attributes;
    return 1;
}

sub resolve_plan_attributes {
    my $self = shift;

    my $variant_reporting_plan = $self->plan;
    my $specific_plan = $variant_reporting_plan->get_plan('expert', $self->name);
    for my $name (keys %{$specific_plan->params}) {
        $self->$name($specific_plan->params->{$name});
    }
    my $translations;
    eval { $translations = $self->provider->get_attribute('translations') };
    if (my $error = Exception::Class->caught('NoTranslationsException')) {
        die $error;
    }
    else {
        $self->translate_inputs($translations);
    }
}

sub plan {
    my $self = shift;

    return Genome::VariantReporting::Framework::Plan::MasterPlan->create_from_json($self->plan_json);
}

sub planned_output_names {
    my $self = shift;

    my @properties = $self->__meta__->properties(
        is_output => 1, is_planned => 1);
    return map {$_->property_name} @properties;
}

sub provider {
    my $self = shift;

    return Genome::VariantReporting::Framework::Component::ResourceProvider->create_from_json($self->provider_json);
}

sub provided_output_names {
    my $self = shift;

    my @properties = $self->__meta__->properties(
        is_output => 1, is_provided => 1);
    return map {$_->property_name} @properties;
}

sub resolve_provided_attributes {
    my $self = shift;

    for my $name ($self->provided_output_names) {
        $self->$name($self->provider->get_attribute($name));
    }
}

# TODO this is not covered by tests
sub validate_with_plan_params {
    my ($self, $params) = validate_pos(@_, 1, 1);

    my @errors = $self->__planned_output_errors__($params);
    if (@errors) {
        $self->print_errors(@errors);
        die $self->error_message("Failed to validate_with_plan_params with params:\n" . Data::Dumper::Dumper $params);
    }
    return;
}

sub __planned_output_errors__ {
    my ($self, $params) = validate_pos(@_, 1, 1);
    my $needed = Set::Scalar->new($self->planned_output_names);
    return $self->_get_missing_errors($params, $needed),
        $self->_get_extra_errors($params, $needed);
}

sub __provided_output_errors__ {
    my ($self, $params) = validate_pos(@_, 1, 1);
    my $needed = Set::Scalar->new($self->provided_output_names);
    return $self->_get_missing_errors($params, $needed);
}

sub _get_missing_errors {
    my ($self, $params, $needed) = validate_pos(@_, 1, 1, 1);

    my $have = Set::Scalar->new(keys %{$params});
    my @errors;
    unless($needed->is_equal($have)) {
        if (my $still_needed = $needed - $have) {
            push @errors, UR::Object::Tag->create(
                type => 'error',
                properties => [$still_needed->members],
                desc => sprintf("Parameters required by adaptor (%s) but not provided: (%s)", 
                    $self->class, join(",", $still_needed->members)),
            );
        }
    }

    return @errors;
}

sub _get_extra_errors {
    my ($self, $params, $needed) = validate_pos(@_, 1, 1, 1);

    my $have = Set::Scalar->new(keys %{$params});
    my @errors;
    unless($needed->is_equal($have)) {
        if (my $not_needed = $have - $needed) {
            push @errors, UR::Object::Tag->create(
                type => 'error',
                properties => [$not_needed->members],
                desc => sprintf("Parameters provided but not required by adaptor (%s): (%s)",
                    $self->class, join(",", $not_needed->members)),
            );
        }
    }

    return @errors;
}

1;
