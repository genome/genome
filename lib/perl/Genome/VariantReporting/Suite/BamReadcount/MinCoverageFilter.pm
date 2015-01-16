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
            doc => 'Miminum coverage',
        },
    },
    doc => q{Filter variants that don't meet minimum coverage},
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

    my $readcount_entries = $self->get_readcount_entries($entry);
    unless ($readcount_entries) {
        return $self->pass_all_sample_alts($entry);
    }

    my %return_value;

    for my $alt_allele (@alt_alleles) {
        my $readcount_entry = $readcount_entries->{$alt_allele};
        if (!defined $readcount_entry) {
            $return_value{$alt_allele} = 1;
        }
        elsif ($readcount_entry->depth >= $self->min_coverage) {
            $return_value{$alt_allele} = 1;
        }
        else {
            $return_value{$alt_allele} = 0;
        }
    }
    return %return_value;
}

sub vcf_description {
    my $self = shift;
    return "Coverage for sample " . $self->sample_name . " is larger than or equal to " . $self->min_coverage;
}

sub vcf_id {
    my $self = shift;
    return "MINCOV" . $self->min_coverage . '_' . $self->sample_name;
}

1;
