package Genome::VariantReporting::Generic::NCallersFilter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Generic::NCallersFilter {
    is => ['Genome::VariantReporting::Framework::Component::Filter', 'Genome::VariantReporting::Framework::Component::WithSampleName'],
    has => [
        min_callers => {
            is => 'Integer',
            doc => 'Variant must be called by at least this many callers',
        },
        valid_callers => {
            is => 'String',
            is_many => 1,
            default_value => [qw(VarscanSomatic Sniper Strelka)],
            doc => 'List of variant callers to include in determination for filtering',
        },
    ],
    doc => q{Filter out variants that weren't called by the minimum number of the specified callers},
};

sub name {
    return 'n-callers';
}

sub requires_annotations {
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
        my $sample_name = $self->sample_name_with_suffix($caller_name);
        my $sample_index = eval{ $entry->{header}->index_for_sample_name($sample_name) };
        my $error = $@;
        if ($error =~ /^\QSample name $sample_name not found in header\E/) {
            next;
        }
        my @sample_alt_alleles = $entry->alt_bases_for_sample($sample_index);
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

sub vcf_id {
    my $self = shift;
    return sprintf("MIN_CALLERS_%s_%s", $self->min_callers, $self->sample_name);
}

sub vcf_description {
    my $self = shift;
    return sprintf("Variant was called by at least %s callers for sample %s", $self->min_callers, $self->sample_name);
}

1;
