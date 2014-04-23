package Genome::Annotation::AlleleInGenotypeFilter;

use strict;
use warnings;
use Genome;

class Genome::Annotation::AlleleInGenotypeFilter {
    is => 'Genome::Annotation::FilterBase',
    has => [
        sample_index => {
            is => 'Integer',
        },
    ],
};

sub name {
    return 'allele-in-genotype';
}

sub requires_experts {
    return ();
}

sub process_entry {
    my $self = shift;
    my $entry = shift;

    my %return_values;
    for my $alt_allele (@{$entry->{alternate_alleles}}) {
        $return_values{$alt_allele} = 0;
    }

    my @sample_alt_alleles = sort $entry->alt_bases_for_sample($self->sample_index);
    for my $sample_alt_allele (@sample_alt_alleles) {
        $return_values{$sample_alt_allele} = 1;
    }

    return %return_values;
}

1;

