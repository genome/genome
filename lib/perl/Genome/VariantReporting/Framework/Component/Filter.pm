package Genome::VariantReporting::Framework::Component::Filter;

use strict;
use warnings FATAL => 'all';
use Genome;
use Carp qw(confess);

class Genome::VariantReporting::Framework::Component::Filter {
    is => 'Genome::VariantReporting::Framework::Component::Interpreter',
    is_abstract => 1,
};

sub filter_entry {
    my $self = shift;
    my $class = $self->class;
    confess "Abstract method 'filter_entry' must be defined in class '$class'";
}

sub available_fields {
    return qw/filter_status/;
}

sub interpret_entry {
    my ($self, $entry, $passed_alt_alleles) = @_;

    my %filter_output = $self->filter_entry($entry);

    my %return_values;
    for my $variant_allele (@$passed_alt_alleles) {
        $return_values{$variant_allele} =  {
            'filter_status' => $filter_output{$variant_allele},
        };
    }
    return %return_values;
}

sub available_fields {
    return qw(filter_status);
}


1;
