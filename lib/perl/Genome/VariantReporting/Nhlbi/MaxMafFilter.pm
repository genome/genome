package Genome::VariantReporting::Nhlbi::MaxMafFilter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Nhlbi::MaxMafFilter {
    is => ['Genome::VariantReporting::Nhlbi::ComponentBase', 'Genome::VariantReporting::Framework::Component::Filter'],
    has => [
        max_maf => {
            is => 'Number',
            doc => 'Maximum minor allele frequency',
        },
    ],
};

sub name {
    return 'max-maf'
}

sub requires_annotations {
    return ('nhlbi');
}

sub filter_entry {
    my $self = shift;
    my $entry = shift;
    my %return_values;

    my $maf = $self->get_maf_for_entry($entry);

    if (!defined $maf or $maf <= $self->max_maf) {
        return map {$_=> 1} @{$entry->{alternate_alleles}};
    }
    else {
        return map {$_=> 0} @{$entry->{alternate_alleles}};
    }

    return %return_values;
}

1;

