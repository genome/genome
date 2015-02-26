package Genome::VariantReporting::Framework::Component::Expert::Command;

use strict;
use warnings FATAL => 'all';
use Genome;
use Set::Scalar;
use JSON;
use Params::Validate qw(validate validate_pos :types);
use List::MoreUtils qw(apply);

my $_JSON_CODEC = new JSON->allow_nonref;

use Genome::VariantReporting::Framework::FileLookup qw(
    is_file
    calculate_lookup
);

class Genome::VariantReporting::Framework::Component::Expert::Command {
    is_abstract => 1,
    is => ['Genome::Command::DelegatesToResult', 'Genome::VariantReporting::Framework::Component::WithTranslatedInputs'],
    attributes_have => {
        is_planned => {
            is => "Boolean",
            default => 0,
        },
    },
    has_structural_input => [
        input_vcf => {
            is => 'Path',
        },
        variant_type => {
            is => 'Text',
            valid_values => ['snvs', 'indels'],
            doc => "The type of variant the input_result represents",
        },
        process_id => {
            is => 'Text',
        },
        plan_json => {
            is => 'Text',
        }
    ],
    has_transient_structural_optional => [
        requestor => {
            is => 'Genome::Process',
            id_by => 'process_id',
        },
    ],
    has_optional_structural_output => [
        output_vcf => {
            is => 'Path',
        },
    ],
};

sub name {
    die "Abstract";
}

sub resolve_plan_attributes {
    my $self = shift;

    my $variant_reporting_plan = $self->plan;
    my $specific_plan = $variant_reporting_plan->get_plan('expert', $self->name);
    while (my ($name, $value) = each %{$specific_plan->run_params}) {
        $self->$name($value);
    }
    return;
}

sub plan {
    my $self = shift;

    return Genome::VariantReporting::Framework::Plan::MasterPlan->create_from_json($self->plan_json);
}

sub planned_names {
    my $self = shift;

    my @properties = $self->__meta__->properties(is_planned => 1);
    return map {$_->property_name} @properties;
}

# TODO this is not covered by tests
sub validate_with_plan_params {
    my ($self, $params) = validate_pos(@_, 1, 1);

    my @errors = $self->__planned_errors__($params);
    if (@errors) {
        $self->print_errors(@errors);
        die $self->error_message("Failed to validate_with_plan_params with params:\n" . Data::Dumper::Dumper $params);
    }
    return;
}

sub __planned_errors__ {
    my ($self, $params) = validate_pos(@_, 1, 1);
    my $needed = Set::Scalar->new($self->planned_names);
    return Genome::VariantReporting::Framework::Utility::get_missing_errors($self->class, $params, $needed, "Parameters", "run"),
        $self->_get_extra_errors($params, $needed);
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
                desc => sprintf("Parameters provided but not required by expert (%s): (%s)",
                    $self->class, join(",", $not_needed->members)),
            );
        }
    }

    return @errors;
}

sub result_class {
    my $self = shift;
    my $class = $self->class;
    die "Abstract method 'result_class' must be defined in class $class";
}

sub post_get_or_create {
    my $self = shift;
    $self->output_vcf($self->output_result->output_file_path);
    return 1;
}

sub input_names {
    my $self = shift;
    return ($self->is_many_input_names, $self->is_not_many_input_names);
}

sub is_many_input_names {
    my $self = shift;

    return apply {s/_lookup$//} $self->result_class->is_many_property_names;
}

sub is_not_many_input_names {
    my $self = shift;

    return apply {s/_lookup$//} $self->result_class->is_not_many_property_names;
}

sub input_hash {
    my $self = shift;

    $self->resolve_plan_attributes;

    my %hash;
    for my $input_name ($self->is_many_input_names) {
        next unless $self->can($input_name);
        my $value = [$self->$input_name];
        $hash{$input_name} = $value;
        if (is_file($value->[0])) {
            $hash{$input_name . '_lookup'} = [map {calculate_lookup($_)} @{$value}];
        }
    }
    for my $input_name ($self->is_not_many_input_names) {
        next unless $self->can($input_name);
        my $value = $self->$input_name;
        if (is_hashref($value)) {
            $hash{$input_name} = json_encode($value);
        } else {
            $hash{$input_name} = $value;
        }

        if (is_file($value)) {
            $hash{$input_name . '_lookup'} = calculate_lookup($self->$input_name);
        }
    }

    $hash{test_name} = $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME};

    return %hash;
}

sub is_hashref {
    my $value = shift;

    if (ref $value eq 'HASH') {
        return 1;
    } else {
        return 0;
    }
}

sub json_encode {
    my $value = shift;

    return $_JSON_CODEC->canonical->encode($value);
}
