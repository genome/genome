package Genome::VariantReporting::Suite::Vep::Tier1Filter;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::VepConsequenceParser;

class Genome::VariantReporting::Suite::Vep::Tier1Filter {
    is => ['Genome::VariantReporting::Framework::Component::Filter'],
    has => [
    ],
    doc => q{Filter out variants that aren't considered tier 1},
};

sub name {
    return 'tier1';
}

sub requires_annotations {
    return ('vep');
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
        transcript_ablation
        splice_donor_variant
        splice_acceptor_variant
        stop_gained
        frameshift_variant
        stop_lost
        initiator_codon_variant
        inframe_insertion
        inframe_deletion
        missense_variant
        transcript_amplification
        splice_region_variant
        incomplete_terminal_codon_variant
        synonymous_variant
        stop_retained_variant
        coding_sequence_variant
        mature_miRNA_variant
    );

    my %consequences;
    map { $consequences{$_} = 1 } @types;
    return \%consequences;
}

sub vcf_id {
    return 'TIER1';
}

sub vcf_description {
    return 'Variant is a tier 1 variant type';
}

1;
