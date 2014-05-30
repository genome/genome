package Genome::VariantReporting::Filter::NCallersFilter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Filter::NCallersFilter {
    is => ['Genome::VariantReporting::Filter::Base', 'Genome::VariantReporting::Filter::WithSampleName'],
    has => [
        min_callers => {
            is => 'Integer',
            doc => 'Variant must be called by at least this many callers',
        },
        valid_callers => {
            is => 'String',
            is_many => 1,
            default_value => [qw(VarscanSomatic Sniper Strelka)],
        },
    ],
};

sub name {
    return 'n-callers';
}

sub requires_experts {
    return ();
}

sub filter_entry {
    my $self = shift;
    my $entry = shift;

    my %return_values;
    my %caller_counts;
    for my $alt_allele (@{$entry->{alternate_alleles}}) {
        $return_values{$alt_allele} = 0;
        $caller_counts{$alt_allele} = 0;
    }

    for my $caller_name ($self->valid_callers) {
        my $sample_index = $entry->{header}->index_for_sample_name(
                                      $self->sample_name_with_suffix($caller_name));
        my @sample_alt_alleles = eval{ $entry->alt_bases_for_sample($sample_index)};
        for my $sample_alt_allele (@sample_alt_alleles) {
            $caller_counts{$sample_alt_allele}++;
        }
    }
    for my $alt_allele (keys %return_values) {
        if ($caller_counts{$alt_allele} >= $self->min_callers) {
            $return_values{$alt_allele} = 1;
        }
    }

    return %return_values;
}
1;

