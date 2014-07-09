package Genome::VariantReporting::Framework::Test::Interpreter;

use strict;
use warnings FATAL => 'all';
use Genome;
use Data::Dump qw(pp);

class Genome::VariantReporting::Framework::Test::Interpreter {
    is => 'Genome::VariantReporting::Framework::Component::Interpreter',
};

sub name {
    return '__test__';
}

sub requires_annotations {
    return qw(__test__);
}

sub interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;
    for my $variant_allele (@$passed_alt_alleles) {
        $return_values{$variant_allele}->{info} = pp($entry->info_for_allele($variant_allele));
    }

    return %return_values;
}

1;
