package Genome::Annotation::PositionInterpreter;

use strict;
use warnings;
use Genome;

class Genome::Annotation::PositionInterpreter {
    is => 'Genome::Annotation::InterpreterBase',
};

sub name {
    return 'position';
}

sub process_entry {
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
