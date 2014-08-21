package Genome::VariantReporting::Generic::HomopolymerFilter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Generic::HomopolymerFilter {
    is  => ['Genome::VariantReporting::Framework::Component::Filter'],
    has => [
        info_tag => {
            is  => 'String',
            doc => 'custom tag name in the info field to show homopolymer status, like HOMO=0,0,1',
        },
        ],
};

sub name {
    return 'homopolymer';
}

sub requires_annotations {
    return ();
}

sub filter_entry {
    my ($self, $entry) = @_;
    my $tag = $self->info_tag;

    unless (exists $entry->{header}->info_types->{$tag}) {
        die sprintf("VCF header does not contain info tag (%s)", $tag);
    }

    my %return_values;

    for my $alt_allele (@{$entry->{alternate_alleles}}) {
        $return_values{$alt_allele} = $entry->info_for_allele($alt_allele, $tag);
    }

    return %return_values;
}


1;
