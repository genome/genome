package Genome::VariantReporting::ComponentBase;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::ComponentBase {
    is_abstract => 1,
};

sub validate {
    my $self = shift;

    my @errors = $self->__errors__;
    if (@errors) {
        $self->print_errors(@errors);
        die $self->error_message("Failed to validate");
    }
    return;
}

sub print_errors {
    my ($self, @errors) = @_;

    for my $error (@errors) {
        my @properties = $error->properties;
        $self->error_message("Property " .
            join(',', map { "'$_'" } @properties) .
            ': ' . $error->desc);
    }
    return;
}

sub part {
    my $self = shift;
    return (split(/::/, $self->class))[-1];
}

1;
