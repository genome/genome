package Genome::VariantReporting::Generic::RsidInterpreter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Generic::RsidInterpreter {
    is => 'Genome::VariantReporting::Framework::Component::Interpreter',
    doc => 'Output a list of rsids',
};

sub name {
    return 'rsid';
}

sub requires_annotations {
    ();
}

sub field_descriptions {
    return (
        rsid => "A list of rsids at this position",
    );
}

sub _interpret_entry {
    my $self = shift;
    my $entry = shift;

    my $passed_alt_alleles = shift;

    my %return_values;
    for my $variant_allele (@$passed_alt_alleles) {
        $return_values{$variant_allele} = {
            rsid => join(",", @{$entry->{identifiers}}),
        };
    }
    return %return_values;
}

1;

