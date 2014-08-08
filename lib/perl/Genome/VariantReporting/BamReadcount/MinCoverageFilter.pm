package Genome::VariantReporting::BamReadcount::MinCoverageFilter;

use strict;
use warnings FATAL => 'all';
use Genome;
use List::Util qw/first/;

class Genome::VariantReporting::BamReadcount::MinCoverageFilter {
    is => ['Genome::VariantReporting::Framework::Component::Filter', 'Genome::VariantReporting::BamReadcount::ComponentBase'],
    has => {
        min_coverage => {
            is => 'Number',
        },
    },
};

sub name {
    return 'min-coverage';
}

sub requires_annotations {
    return ('bam-readcount');
}

sub filter_entry {
    my ($self, $entry) = @_;

    my @alt_alleles = @{$entry->{alternate_alleles}};

    my $readcount_entry = $self->get_readcount_entry($entry);
    unless ($readcount_entry) {
        return $self->pass_all_sample_alts($entry);
    }

    if ($readcount_entry->depth >= $self->min_coverage) {
        return map { $_ => 1 } @alt_alleles;
    }
    return map { $_ => 0 } @alt_alleles;
}

1;
