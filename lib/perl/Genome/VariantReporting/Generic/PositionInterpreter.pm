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

sub _interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;
    for my $variant_allele (@$passed_alt_alleles) {
        my ($stop, $reference, $variant);
        if (length($variant_allele) < length($entry->{reference_allele})) {
            ($stop, $reference, $variant) = _interpret_deletion($entry, $variant_allele);
        }
        elsif (length($variant_allele) > length($entry->{reference_allele})) {
            ($stop, $reference, $variant) = _interpret_insertion($entry, $variant_allele);
        }
        else {
            $stop = $entry->{position} + length($variant_allele) - 1;
            $reference = $entry->{reference_allele};
            $variant = $variant_allele;
        }
        $return_values{$variant_allele} = {
            chromosome_name => $entry->{chrom},
            start           => $entry->{position},
            stop            => $stop,
            reference       => $reference,
            variant         => $variant,
        };
    }

    return %return_values;
}

sub _interpret_deletion {
    my ($entry, $variant_allele) = @_;
    my $deletion_length = length($entry->{reference_allele}) - length($variant_allele);
    my $stop = $entry->{position} + $deletion_length;
    my $reference = substr($entry->{reference_allele}, length($variant_allele));
    my $variant = "-";
    return ($stop, $reference, $variant);
}

sub _interpret_insertion {
    my ($entry, $variant_allele) = @_;
    my $stop = $entry->{position} + 1;
    my $reference = "-";
    my $variant = substr($variant_allele, length($entry->{reference_allele}));
    return ($stop, $reference, $variant);
}

1;
