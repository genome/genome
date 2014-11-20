package Genome::VariantReporting::Generic::ContainsTagFilter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Generic::ContainsTagFilter {
    is => 'Genome::VariantReporting::Framework::Component::Filter',
    has => [
        info_tag => {
            is => "String",
            doc => "Entry must contain this tag in the info field to keep",
        },
    ],
    doc => q{Filter out variants that don't contain the specified info tag},
};

sub name {
    return 'contains-tag';
}

sub requires_annotations {
    return ();
}

sub filter_entry {
    my $self = shift;
    my $entry = shift;

    unless (exists $entry->{header}->info_types->{$self->info_tag}) {
        die sprintf("VCF header does not contain info tag (%s)", $self->info_tag);
    }

    my %return_values;

    my $contains_tag;
    if (defined $entry->info($self->info_tag)) {
        $contains_tag = 1;
    } else {
        $contains_tag = 0;
    }

    for my $alt_allele (@{$entry->{alternate_alleles}}) {
        $return_values{$alt_allele} = $contains_tag;
    }

    return %return_values;
}

sub vcf_id {
    my $self = shift;
    return 'CONTAINS_TAG_' . $self->info_tag;
}

sub vcf_description {
    my $self = shift;
    return 'INFO field for entry contains tag ' . $self->info_tag;
}

1;

