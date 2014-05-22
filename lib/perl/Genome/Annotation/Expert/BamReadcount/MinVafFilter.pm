package Genome::Annotation::Expert::BamReadcount::MinVafFilter;

use strict;
use warnings;
use Genome;
use Genome::Annotation::Expert::BamReadcount::VafCalculator;

class Genome::Annotation::Expert::BamReadcount::MinVafFilter {
    is => ['Genome::Annotation::Filter::Base', 'Genome::Annotation::Expert::BamReadcount::ComponentBase'],
    has => [
        min_vaf => {
            is => 'Number',
        },
    ],
};

sub name {
    return 'min-vaf';
}

sub requires_experts {
    return ('bam-readcount');
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
        my $vaf = Genome::Annotation::Expert::BamReadcount::VafCalculator::calculate_vaf(
            $self->get_readcount_entry($entry), $sample_alt_allele);
        if ($vaf >= $self->min_vaf) {
            $return_values{$sample_alt_allele} = 1;
        }
    }

    return %return_values;
}

1;

