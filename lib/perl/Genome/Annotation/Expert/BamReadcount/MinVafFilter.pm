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
    my %vafs = Genome::Annotation::Expert::BamReadcount::VafCalculator::calculate_vaf_for_multiple_alleles(
        $self->get_readcount_entry($entry), \@sample_alt_alleles);
    for my $allele (keys %vafs) {
        if ($vafs{$allele} >= $self->min_vaf) {
            $return_values{$allele} = 1;
        }
    }

    return %return_values;
}

1;

