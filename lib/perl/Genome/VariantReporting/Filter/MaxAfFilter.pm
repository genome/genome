package Genome::VariantReporting::Filter::MaxAfFilter;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::DbsnpAFParser;

class Genome::VariantReporting::Filter::MaxAfFilter {
    is => 'Genome::VariantReporting::Component::Filter',
    has => [
        max_af => {
            is => 'Number',
            doc => 'Maximum allele frequency',
        },
    ],
};

sub name {
    return 'max-af'
}

sub requires_experts {
    return ('dbsnp');
}

sub filter_entry {
    my $self = shift;
    my $entry = shift;
    my %return_values;

    my $parser = _caf_parser($entry->{header});
    my $caf = $parser->process_entry($entry);

    if (!defined $caf) {
        return map {$_ => 1} @{$entry->{alternate_alleles}};
    }

    for my $alt_allele (@{$entry->{alternate_alleles}}) {
        if ($caf->{$alt_allele} <= $self->max_af) {
            $return_values{$alt_allele} = 1;
        }
        else {
            $return_values{$alt_allele} = 0;
        }
    }

    return %return_values;
}

sub _caf_parser {
    my $header = shift;
    return Genome::File::Vcf::DbsnpAFParser->new($header);
}

Memoize::memoize('_caf_parser');

1;

