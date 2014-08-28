package Genome::VariantReporting::Joinx::ThousandGenomes::MaxAfFilter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Joinx::ThousandGenomes::MaxAfFilter {
    is => ['Genome::VariantReporting::Framework::Component::Filter'],
    has => [
        max_af => {
            is => 'Number',
            doc => 'Maximum minor allele frequency',
        },
    ],
};

sub name {
    return '1kg-max-af'
}

sub requires_annotations {
    return ('1kg');
}

sub filter_entry {
    my $self = shift;
    my $entry = shift;
    my %return_values;

    my $af = $entry->info("AF");

    if (!defined $af or $af <= $self->max_af) {
        return map {$_=> 1} @{$entry->{alternate_alleles}};
    }
    else {
        return map {$_=> 0} @{$entry->{alternate_alleles}};
    }

    return %return_values;
}

1;

