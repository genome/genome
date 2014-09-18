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

sub field_descriptions {
    return (
        chromosome_name => 'Chromosome: an identifier from the reference genome or an angle-bracketed ID String ("<ID>") pointing to a contig in the assembly file',
        start => 'The start position of the variant. One-based',
        stop => 'The end position of the variant. One-based, inclusive',
        reference => 'The reference allele at this position',
        variant => 'The variant called at this position',
    );
}

sub _interpret_entry {
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
