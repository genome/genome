package Genome::Annotation::Reporter::WithHeader;

use strict;
use warnings FATAL => 'all';
use Genome;
use Memoize qw();

class Genome::Annotation::Reporter::WithHeader {
    is => 'Genome::Annotation::Reporter::Base',
    is_abstract => 1,
    has => {
        null_character => {
            is => 'Text',
            default => '-'
        },
        delimiter => {
            is => 'Text',
            default => "\t",
        },
    }
};

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__;
    my %available_fields = $self->available_fields_dict;
    for my $header ($self->headers) {
        unless(defined($available_fields{$header})) {
                die $self->error_message("Interpreter field for $header is not defined. Do you need to overwrite available_fields_dict to provide the correct mapping?");
        }
    }
}

sub requires_interpreters_classes {
    my $self = shift;
    my @interpreters;
    for my $interpreter_name ($self->requires_interpreters) {
        push @interpreters, Genome::Annotation::Factory->_load('interpreters')->{$interpreter_name};
    }
    return @interpreters;
}

sub initialize {
    my $self = shift;
    my $output_dir = shift;

    $self->SUPER::initialize($output_dir);
    $self->print_headers();
}

sub headers {
    die "abstract";
}

sub print_headers {
    my $self = shift;

    my @headers = $self->headers();
    $self->_output_fh->print(join($self->delimiter, @headers) . "\n");
}

# Default dictionary that maps headers to interpreter fields
# This can be used if the headers are named exactly the same as the
# interpreters' fields. The interpreter fields need to be unique as well.
# This also works if you are only interested in a subset of the interpreters'
# fields as long as the above requirements are met.
# Overwrite in child class if different behavior desired.
sub available_fields_dict {
    my $self = shift;

    my @interpreters = $self->requires_interpreters_classes;
    my %available_fields;
    for my $interpreter (@interpreters) {
        for my $field ($interpreter->available_fields()) {
            if (defined $available_fields{$field}) {
                die $self->error_message("Fields are not unique. Field: %s, Interpreters: %s and %s",
                    $field, $interpreter->name, $available_fields{$field}->{interpreter});
            }
            $available_fields{$field} = {
                interpreter => $interpreter->name,
                field => $field,
            }
        }
    }
    return %available_fields;
}
Memoize::memoize('available_fields_dict');

# Default report method
# Prints the fields in order of the headers.
# Overwrite in child class if different behavior desired.
sub report {
    my $self = shift;
    my $interpretations = shift;

    my %fields = $self->available_fields_dict();
    for my $allele (keys %{$interpretations->{($self->requires_interpreters)[0]}}) {
        for my $header ($self->headers()) {
            my $interpreter = $fields{$header}->{interpreter};
            my $field = $fields{$header}->{field};
            $self->_output_fh->print($self->_format($interpretations->{$interpreter}->{$allele}->{$field}) . "\t");
        }
        $self->_output_fh->print("\n");
    }
}

sub _format {
    my $self = shift;
    my $string = shift;

    if (defined $string ) {
        return $string;
    }
    else {
        return $self->null_character;
    }
}

1;
