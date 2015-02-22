package Genome::VariantReporting::Suite::Joinx::Homopolymer::HomopolymerFilter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Suite::Joinx::Homopolymer::HomopolymerFilter {
    is  => ['Genome::VariantReporting::Framework::Component::Filter'],
    has => [
        info_tag => {
            is  => 'String',
            doc => 'custom tag name in the info field to show homopolymer status, like HOMP_FILTER=0,1',
        },
    ],
    doc => 'Filter out indels hitting homopolymers',
};

sub name {
    return 'homopolymer';
}

sub requires_annotations {
    return ('homopolymer');
}

sub filter_entry {
    my ($self, $entry) = @_;
    my $tag = $self->info_tag;

    unless (exists $entry->{header}->info_types->{$tag}) {
        die sprintf("VCF header does not contain info tag (%s)", $tag);
    }

    my %return_values;

    for my $alt_allele (@{$entry->{alternate_alleles}}) {
        my $status = $entry->info_for_allele($alt_allele, $tag);
        if (defined $status) {
            $status =~ tr/01/10/;   #joinx sets 1 for homopolymer hit, 0 for not hit
            $return_values{$alt_allele} = $status;
        }
        else { #No HOMP_FILTER in INFO column
            $return_values{$alt_allele} = 1;
        }
    }

    return %return_values;
}

sub vcf_id {
    return 'HOMOPOLYMERS';
}

sub vcf_description {
    return 'Filter out indels hitting homopolymers';
}


1;
