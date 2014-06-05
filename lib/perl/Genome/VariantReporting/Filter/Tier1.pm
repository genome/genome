package Genome::VariantReporting::Filter::Tier1;

use strict;
use warnings;
use Genome;
use Scalar::Util qw(looks_like_number);
use Genome::File::Vcf::VepConsequenceParser;

class Genome::VariantReporting::Filter::Tier1 {
    is => ['Genome::VariantReporting::Filter::Base'],
    has => [
    ],
};

sub name {
    return 'tier1';
}

sub filter_entry {
    my $self = shift;
    my $entry = shift;

    my %return_values;
    my $vep_parser = new Genome::File::Vcf::VepConsequenceParser($entry->{header});
    my $tier1_consequence_types = $self->tier1_consequence_types;

    for my $alt_allele ( @{$entry->{alternate_alleles}} ) {
        my ($transcript) = $vep_parser->transcripts($entry, $alt_allele);
        my $consequence = $transcript->{'consequence'};
        if (defined $consequence and exists $tier1_consequence_types->{lc($consequence)}) {
            $return_values{$alt_allele} = 1;
        } else {
            $return_values{$alt_allele} = 0;
        }
    }

    return %return_values;
}

sub tier1_consequence_types {
    my @types = qw(
        synonymous_variant
        missense_variant
        inframe_insertion
        inframe_deletion
        stop_gained
        frameshift_variant
        coding_sequence_variant
        stop_lost
        stop_retained_variant
        incomplete_terminal_codon_variant
    );

    my %consequences;
    map { $consequences{$_} = 1 } @types;
    return \%consequences;
}

1;
