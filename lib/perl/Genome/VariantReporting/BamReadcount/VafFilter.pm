package Genome::VariantReporting::BamReadcount::VafFilter;

use strict;
use warnings;
use Genome;
use Genome::VariantReporting::BamReadcount::VafCalculator;

class Genome::VariantReporting::BamReadcount::VafFilter {
    is => ['Genome::VariantReporting::Framework::Component::Filter', 'Genome::VariantReporting::BamReadcount::ComponentBase'],
    has_optional => [
        min_vaf => {
            is => 'Number',
            doc => 'The inclusive lower bound for VAF values that will pass the filter (will be kept)',
        },
        max_vaf => {
            is => 'Number',
            doc => 'The inclusive upper bound for VAF values that will pass the filter (will be kept)',
        },
    ],
};

sub name {
    return 'vaf-cutoff';
}

sub requires_experts {
    return ('bam-readcount');
}

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;

    unless (defined($self->min_vaf) || defined($self->max_vaf)) {
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => [$self->min_vaf, $self->max_vaf],
            desc => "Must define at least one of min_vaf or max_vaf",
        );
    }

    if (defined($self->min_vaf) && defined($self->max_vaf)
        && $self->min_vaf > $self->max_vaf
    ) {
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => [$self->min_vaf, $self->max_vaf],
            desc => "Max_vaf must be larger or equal to min_vaf",
        );

    }

    return @errors;
}

sub filter_entry {
    my $self = shift;
    my $entry = shift;
    my %return_values;
    for my $alt_allele (@{$entry->{alternate_alleles}}) {
        $return_values{$alt_allele} = 0;
    }
    my %vafs = Genome::VariantReporting::BamReadcount::VafCalculator::calculate_vaf_for_all_alts(
        $entry,
        $self->get_readcount_entry($entry));
    for my $allele (keys %vafs) {
        if ($self->passes_filter($vafs{$allele})) {
            $return_values{$allele} = 1;
        }
    }

    return %return_values;
}

sub passes_filter {
    my $self = shift;
    my $vaf_value = shift;

    if (defined($self->min_vaf)) {
        unless ($vaf_value >= $self->min_vaf) {
            return 0;
        }
    }

    if (defined($self->max_vaf)) {
        unless ($vaf_value <= $self->max_vaf) {
            return 0;
        }
    }

    return 1;
}

1;

