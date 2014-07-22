package Genome::VariantReporting::Generic::PositionInterpreter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Generic::PositionInterpreter {
    is => 'Genome::VariantReporting::Framework::Component::Interpreter',
};

sub name {
    return 'position';
}

sub requires_annotations {
    ();
}

sub available_fields {
    return qw/
        chromosome_name
        start
        stop
        reference
        variant
    /;
}

sub process_interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;
    for my $variant_allele (@$passed_alt_alleles) {
        $return_values{$variant_allele} = {
            chromosome_name => $entry->{chrom},
            start           => $entry->{position},
            stop            => $entry->{position} + length($variant_allele) - 1,
            reference       => $entry->{reference_allele},
            variant         => $variant_allele
        };
    }

    return %return_values;
}

1;
