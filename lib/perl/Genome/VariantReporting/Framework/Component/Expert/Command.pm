package Genome::VariantReporting::Framework::Component::Expert::Command;

use strict;
use warnings FATAL => 'all';
use Genome;
use Set::Scalar;
use JSON;

my $_JSON_CODEC = new JSON->allow_nonref;

use Genome::VariantReporting::Framework::FileLookup qw(
    is_file
    calculate_lookup
);

class Genome::VariantReporting::Framework::Component::Expert::Command {
    is_abstract => 1,
    is => ['Genome::Command::DelegatesToResult', 'Genome::VariantReporting::Framework::Component::Base'],
    has_input => [
        input_vcf => {
            is => 'Path',
        },
        variant_type => {
            is => 'Text',
            valid_values => ['snvs', 'indels'],
            doc => "The type of variant the input_result represents",
        },
        user => {
            is => 'Genome::Process',
            is_optional => 1,
            id_by => 'process_id',
        },
        process_id => {
            is => 'Text',
        },
    ],
    has_optional_output => [
        output_vcf => {
            is => 'Path',
        },
        output_result => {
            is => 'Genome::SoftwareResult',
        },
    ],
};

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

    my @properties = $self->__meta__->properties(
        is_many => 1, is_input => 1);
    return map {$_->property_name} @properties;
}

sub is_not_many_input_names {
    my $self = shift;

    my @properties = $self->__meta__->properties(
        is_many => 0, is_input => 1);
    return map {$_->property_name} @properties;
}

sub input_hash {
    my $self = shift;

    my %hash;
    for my $input_name ($self->is_many_input_names) {
        my $value = [$self->$input_name];
        $hash{$input_name} = $value;
        if (is_file($value->[0])) {
            $hash{$input_name . '_lookup'} = [map {calculate_lookup($_)} @{$value}];
        }
    }
    for my $input_name ($self->is_not_many_input_names) {
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
    delete $hash{user};
    delete $hash{process_id};
    delete $hash{label};
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
