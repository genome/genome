package Genome::Annotation::ReporterBase;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::ReporterBase {
    is => 'Genome::Annotation::ComponentBase',
    is_abstract => 1,
};

sub name {
    die "abstract";
}

sub validate {
    my $self = shift;
    if (my @errors = $self->__errors__) {
        $self->print_errors(@errors);
        die $self->error_message("Report (%s) failed validation", $self->name);
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
