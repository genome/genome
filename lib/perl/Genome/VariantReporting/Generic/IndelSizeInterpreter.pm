package Genome::VariantReporting::Generic::IndelSizeInterpreter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Generic::IndelSizeInterpreter {
    is => 'Genome::VariantReporting::Framework::Component::Interpreter',
};

sub name {
    return 'indel-size';
}

sub requires_annotations {
    return qw/ /;
}

sub field_descriptions {
    return (
        indel_size => 'Length of indel',
    );
}

sub _interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;
    for my $alt_allele ( @$passed_alt_alleles ) {
        my $indel_size = abs( length($entry->{reference_allele}) - length($alt_allele) );
        $return_values{$alt_allele} = {
            indel_size => $indel_size
        };
    }

    return %return_values;
}

1;

