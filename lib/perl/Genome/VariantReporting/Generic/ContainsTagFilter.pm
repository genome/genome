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
};

sub name {
    return 'contains-tag';
}

sub requires_experts {
    return ();
}

sub available_fields {
    return qw/contains_tag/;
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

sub interpret_entry {
    my ($self, $entry, $passed_alt_alleles) = @_;

    my %filter_output = $self->filter_entry($entry);

    my %return_values;
    for my $variant_allele (@$passed_alt_alleles) {
        $return_values{$variant_allele} =  {
            contains_tag => $filter_output{$variant_allele},
        };
    }
    return %return_values;
}

1;

