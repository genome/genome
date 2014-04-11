package Genome::Annotation::Detail::Command;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Detail::Command {
    is => 'Command::V2',
};

sub input_names {
    my $self = shift;

    my @properties = $self->__meta__->properties(class_name => $self->class);
    return map {$_->property_name} grep {$_->is_input} @properties;
}

sub input_hash {
    my $self = shift;

    my %hash;
    for my $input_name ($self->input_names) {
        $hash{$input_name} = $self->$input_name;
    }
    return %hash;
}
