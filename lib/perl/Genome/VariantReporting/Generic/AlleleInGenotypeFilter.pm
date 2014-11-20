package Genome::VariantReporting::Generic::AlleleInGenotypeFilter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Generic::AlleleInGenotypeFilter {
    is => ['Genome::VariantReporting::Framework::Component::Filter', 'Genome::VariantReporting::Framework::Component::WithSampleName'],
    has => [
    ],
    doc => q{Filter out alternate alleles that aren't part of the specified sample's genotype},
};

sub name {
    return 'allele-in-genotype';
}

sub requires_annotations {
    return ();
}

sub filter_entry {
    my $self = shift;
    my $entry = shift;

    my %return_values;
    for my $alt_allele (@{$entry->{alternate_alleles}}) {
        $return_values{$alt_allele} = 0;
    }

    my @sample_alt_alleles = sort $entry->alt_bases_for_sample($self->sample_index($entry->{header}));
    for my $sample_alt_allele (@sample_alt_alleles) {
        $return_values{$sample_alt_allele} = 1;
    }

    return %return_values;
}

sub vcf_id {
    my $self = shift;
    return 'ALLELE_IN_SAMPLE_GENOTYPE_' . $self->sample_name;
}

sub vcf_description {
    my $self = shift;
    return 'Filter out alleles that are not part of the genotype of sample ' . $self->sample_name;
}

1;

