package Genome::Annotation::ComponentBase;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::ComponentBase {
    is_abstract => 1,
};

sub validate {
    my $self = shift;
    $self->validate_params;
}

sub validate_params {
    my $self = shift;

    if (my @errors = $self->__errors__) {
        $self->print_errors(@errors);
        die $self->error_message("%s (%s) failed validation", $self->part, $self->name);
    }
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


1;
