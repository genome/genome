package Genome::VariantReporting::Suite::Joinx::ThousandGenomes::MaxAfFilter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Suite::Joinx::ThousandGenomes::MaxAfFilter {
    is => ['Genome::VariantReporting::Framework::Component::Filter'],
    has => [
        max_af => {
            is => 'Number',
            doc => 'Maximum minor allele frequency',
        },
    ],
    doc => 'Filter variants that exceed maximum minor allele frequency',
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

sub vcf_id {
    my $self = shift;
    return 'MAXAF' . $self->max_af;
}

sub vcf_description {
    my $self = shift;
    return 'Filter out variants with minor allele frequency being greater than ' . $self->max_af;
}

1;

