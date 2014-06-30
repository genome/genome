package Genome::VariantReporting::Generic::RsidInterpreter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Generic::RsidInterpreter {
    is => 'Genome::VariantReporting::Framework::Component::Interpreter',
};

sub name {
    return 'rsid';
}

sub requires_experts {
    ();
}

sub available_fields {
    return qw/
        rsid
    /;
}

sub interpret_entry {
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

