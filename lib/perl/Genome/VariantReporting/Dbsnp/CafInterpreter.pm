package Genome::VariantReporting::Dbsnp::CafInterpreter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Dbsnp::CafInterpreter {
    is => ['Genome::VariantReporting::Dbsnp::ComponentBase', 'Genome::VariantReporting::Framework::Component::Interpreter'],
};

sub name {
    return 'caf'
}

sub requires_annotations {
    return qw/
        dbsnp
    /;
}

sub available_fields {
    return qw/
        caf
    /;
}

sub _interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;

    for my $variant_allele (@$passed_alt_alleles) {
        if (!defined $entry->info("CAF")) {
            $return_values{$variant_allele}->{caf} = undef;
        }
        else {
            $return_values{$variant_allele}->{caf} = $entry->info("CAF");
        }
    }

    return %return_values;
}

1;

