package Genome::Annotation::Filter::Base;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Filter::Base {
    is => 'Genome::Annotation::Interpreter::Base',
    is_abstract => 1,
};

sub filter_entry {
    my $self = shift;
    my $class = $self->class;
    die "Abstract method 'filter_entry' must be defined in class '$class'";
}

sub interpret_entry {
    my ($self, $entry, $passed_alt_alleles) = shift;

    my %filter_output = $self->filter_entry($entry);

    my %return_values;
    for my $variant_allele (@$passed_alt_alleles) {
        $return_values{$variant_allele} =  {
            'filter_status' => $filter_output{$variant_allele},
        };
    }
    return %return_values;
}


1;
