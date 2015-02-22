package Genome::VariantReporting::Suite::Joinx::Dbsnp::MaxAfFilter;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::DbsnpAFParser;

class Genome::VariantReporting::Suite::Joinx::Dbsnp::MaxAfFilter {
    is => ['Genome::VariantReporting::Suite::Joinx::Dbsnp::ComponentBase', 'Genome::VariantReporting::Framework::Component::Filter'],
    has => [
        max_af => {
            is => 'Number',
            doc => 'Maximum allele frequency',
        },
    ],
    doc => 'Filter variants that exceed maximum allele frequency cutoff',
};

sub name {
    return 'max-af'
}

sub requires_annotations {
    return ('dbsnp');
}

sub filter_entry {
    my $self = shift;
    my $entry = shift;
    my %return_values;

    my $parser = $self->_caf_parser($entry->{header});
    my $caf = $parser->process_entry($entry);

    if (!defined $caf) {
        return map {$_ => 1} @{$entry->{alternate_alleles}};
    }

    for my $alt_allele (@{$entry->{alternate_alleles}}) {
        if (!defined $caf->{$alt_allele}
            or $caf->{$alt_allele} eq '.'
            or $caf->{$alt_allele} <= $self->max_af) {
            $return_values{$alt_allele} = 1;
        }
        else {
            $return_values{$alt_allele} = 0;
        }
    }

    return %return_values;
}

sub vcf_id {
    my $num = shift->max_af;
    $num =~ s/\.//g;
    return 'MAXAF'.$num;
}

sub vcf_description {
    return 'Filter out variants with maf being greater than '.shift->max_af;
}


1;

