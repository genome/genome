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

sub field_descriptions {
    return (
        info => 'Description of info',
    );
}

sub _interpret_entry {
    my ($self, $entry, $passed_alt_alleles) = @_;

    my %return_values;
    for my $variant_allele (@$passed_alt_alleles) {
        $return_values{$variant_allele}->{info} = pp($entry->info_for_allele($variant_allele));
    }

    return %return_values;
}

1;
