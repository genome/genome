package Genome::VariantReporting::Framework::Component::Expert::Command;

use strict;
use warnings FATAL => 'all';
use Genome;
use Data::Dump qw(pp);
use Set::Scalar;
use JSON;

my $_JSON_CODEC = new JSON->allow_nonref;

use Genome::VariantReporting::Framework::FileLookup qw(
    is_file
    calculate_lookup
);

class Genome::VariantReporting::Framework::Component::Expert::Command {
    is_abstract => 1,
    is => ['Command::V2', 'Genome::VariantReporting::Framework::Component::Base'],
    has_input => [
        input_vcf => {
            is => 'Path',
        },
        variant_type => {
            is => 'Text',
            valid_values => ['snvs', 'indels'],
            doc => "The type of variant the input_result represents",
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

sub shortcut {
    my $self = shift;

    $self->debug_message("Attempting to get a %s with arugments %s",
        $self->result_class, pp($self->input_hash));
    my $result = $self->result_class->get_with_lock($self->input_hash);
    if ($result) {
        $self->debug_message("Found existing result (%s) with output_file_path (%s)",
            $result->id, $result->output_file_path);
        $self->output_result($result);
        $self->output_vcf($result->output_file_path);
        return 1;
    } else {
        $self->debug_message("Found no existing result.");
        return 0;
    }
}

sub execute {
    my $self = shift;

    $self->debug_message("Validating inputs");
    $self->validate();

    $self->debug_message("Attempting to get or create a %s with arugments %s",
        $self->result_class, pp({$self->input_hash}));
    my $result = $self->result_class->get_or_create($self->input_hash);
    $self->debug_message("Got or created result (%s) with output file path (%s)",
        $result->id, $result->output_file_path);
    $self->debug_message("Calculated query for result:\n".pp({$result->calculate_query()}));
    $self->output_result($result);
    $self->output_vcf($result->output_file_path);
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
