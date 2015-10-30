package Genome::VariantReporting::Generic::NCallersFilter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Generic::NCallersFilter {
    is => ['Genome::VariantReporting::Framework::Component::Filter', 'Genome::VariantReporting::Generic::NCallersBase'],
    has => [
        min_callers => {
            is => 'Integer',
            doc => 'Variant must be called by at least this many callers',
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
    for my $alt_allele (@{$entry->{alternate_alleles}}) {
        $return_values{$alt_allele} = 0;
    }
    my %callers = $self->get_callers($entry, $entry->{alternate_alleles});
    for my $alt_allele (keys %return_values) {
        if (@{$callers{$alt_allele}} >= $self->min_callers) {
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
