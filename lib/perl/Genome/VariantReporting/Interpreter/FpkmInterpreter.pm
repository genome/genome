package Genome::VariantReporting::Interpreter::FpkmInterpreter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Interpreter::FpkmInterpreter {
    is => ['Genome::VariantReporting::Component::Interpreter', 'Genome::VariantReporting::Component::Interpreter::Fpkm'],
};

sub name {
    return 'fpkm';
}

sub requires_experts {
    ('fpkm');
}

sub available_fields {
    return qw/
        fpkm
    /;
}

sub interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;
    my %fpkm_for_genotype_allele = $self->fpkm_for_genotype_allele($entry);
    for my $variant_allele (@$passed_alt_alleles) {
        $return_values{$variant_allele} = {
            fpkm => $fpkm_for_genotype_allele{$variant_allele},
        };
    }
    return %return_values;
}

1;
