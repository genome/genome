package Genome::Annotation::Detail::Command;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Detail::Command {
    is => 'Command::V2',
};

sub input_names {
    my $self = shift;
    return ($self->is_many_input_names, $self->is_not_many_input_names);
}

sub is_many_input_names {
    my $self = shift;

    my @properties = $self->__meta__->properties(class_name => $self->class,
        is_many => 1, is_input => 1);
    return map {$_->property_name} @properties;
}

sub is_not_many_input_names {
    my $self = shift;

    my @properties = $self->__meta__->properties(class_name => $self->class,
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
