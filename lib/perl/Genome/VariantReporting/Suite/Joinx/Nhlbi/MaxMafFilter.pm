package Genome::VariantReporting::Suite::Joinx::Nhlbi::MaxMafFilter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Suite::Joinx::Nhlbi::MaxMafFilter {
    is => ['Genome::VariantReporting::Suite::Joinx::Nhlbi::ComponentBase', 'Genome::VariantReporting::Framework::Component::Filter'],
    has => [
        max_maf => {
            is => 'Number',
            doc => 'Maximum minor allele frequency',
        },
        population_code => {
            is => 'String',
            valid_values => ['All', 'EU', 'AA'],
            doc => 'Population for which to compare MAF',
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

    my $maf = $self->get_maf_for_entry($entry, $self->population_code);

    if (!defined $maf or $maf <= $self->max_maf) {
        return map {$_=> 1} @{$entry->{alternate_alleles}};
    }
    else {
        return map {$_=> 0} @{$entry->{alternate_alleles}};
    }

    return %return_values;
}

sub vcf_id {
    my $self = shift;
    my $num  = $self->max_maf;
    $num =~ s/\.//g;
    return 'MAXMAF'.uc($self->population_code).$num;
}

sub vcf_description {
    my $self = shift;
    return 'Filter out variants with '.$self->population_code.' maf being greater than '.$self->max_maf;
}


1;

