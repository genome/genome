package Genome::Annotation::Filter::ContainsTagFilter;

use strict;
use warnings;
use Genome;

class Genome::Annotation::Filter::ContainsTagFilter {
    is => 'Genome::Annotation::Filter::Base',
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

sub filter_entry {
    my $self = shift;
    my $entry = shift;

    my %return_values;

    my $contains_tag = defined $entry->info($self->info_tag);

    for my $alt_allele (@{$entry->{alternate_alleles}}) {
        $return_values{$alt_allele} = $contains_tag;
    }

    return %return_values;
}

1;

