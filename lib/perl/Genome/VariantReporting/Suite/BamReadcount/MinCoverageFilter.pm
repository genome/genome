package Genome::VariantReporting::Suite::BamReadcount::MinCoverageFilter;

use strict;
use warnings FATAL => 'all';
use Genome;
use List::Util qw/first/;

class Genome::VariantReporting::Suite::BamReadcount::MinCoverageFilter {
    is => ['Genome::VariantReporting::Framework::Component::Filter', 'Genome::VariantReporting::Suite::BamReadcount::ComponentBase'],
    has => {
        min_coverage => {
            is => 'Number',
        },
    },
    doc => 'Filter variants that don\'t meet minimum coverage',
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

sub vcf_description {
    my $self = shift;
    return "Coverage for sample " . $self->sample_name . " is larger than or equal to " . $self->min_coverage;
}

sub vcf_id {
    my $self = shift;
    return "MINCOV" . $self->min_coverage;
}

1;
