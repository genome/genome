package Genome::VariantReporting::Interpreter::VariantTypeInterpreter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Interpreter::VariantTypeInterpreter {
    is => 'Genome::VariantReporting::Component::Interpreter',
};

sub name {
    return 'variant-type';
}

sub requires_experts {
    ();
}

sub available_fields {
    return qw/
        variant_type
    /;
}

sub interpret_entry {
    my $self = shift;
    my $entry = shift;

    my $passed_alt_alleles = shift;

    my %return_values;
    for my $variant_allele (@$passed_alt_alleles) {
        $return_values{$variant_allele} = {
            variant_type => variant_type($entry, $variant_allele),
        };
    }
    return %return_values;
}

sub variant_type {
    my $entry = shift;
    my $variant_allele = shift;

    if (length($variant_allele) == length($entry->{reference_allele})) {
        if (length($variant_allele) == 1) {
            return 'snp';
        }
        else {
            return 'other';
        }
    }
    else {
        return 'indel';
    }
}

1;

