package Genome::VariantReporting::Expert::CommandBase;

use strict;
use warnings FATAL => 'all';
use Genome;
use Data::Dump qw(pp);

class Genome::VariantReporting::Expert::CommandBase {
    is_abstract => 1,
    is => ['Command::V2', 'Genome::VariantReporting::ComponentBase'],
    has_input => [
        input_result => {
            is => 'Genome::SoftwareResult',
            doc => "The software result created by the previously run command",
        },
        variant_type => {
            is => 'Text',
            valid_values => ['snvs', 'indels'],
            doc => "The type of variant the input_result represents",
        },
    ],
    has_optional_output => [
        output_result => {
            is => 'Genome::VariantReporting::Expert::ResultBase',
            doc => 'The software result created during command execution',
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
        $self->output_result($result);
        return 1;
    } else {
        return 0;
    }
}

sub execute {
    my $self = shift;

    $self->debug_message("Validating inputs");
    $self->validate();

    $self->debug_message("Attempting to get or create a %s with arugments %s",
        $self->result_class, pp($self->input_hash));
    $self->output_result($self->result_class->get_or_create($self->input_hash));
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
        $hash{$input_name} = [$self->$input_name];
    }
    for my $input_name ($self->is_not_many_input_names) {
        $hash{$input_name} = $self->$input_name;
    }
    return %hash;
}
