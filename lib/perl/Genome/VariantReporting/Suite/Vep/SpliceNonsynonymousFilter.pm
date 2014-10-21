package Genome::VariantReporting::Suite::Vep::SpliceNonsynonymousFilter;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::VepConsequenceParser;
use Genome::VariantReporting::Suite::Vep::SpliceNonsynonymousList;

class Genome::VariantReporting::Suite::Vep::SpliceNonsynonymousFilter {
    is  => ['Genome::VariantReporting::Framework::Component::Filter'],
};


sub name {
    return 'splice-nonsynonymous';
}

sub requires_annotations {
    return ('vep');
}


sub filter_entry {
    my ($self, $entry) = @_;
    my %return_values;

    my $vep_parser = Genome::File::Vcf::VepConsequenceParser->new($entry->{header});

    for my $alt_allele (@{$entry->{alternate_alleles}}) {
        my ($transcript) = $vep_parser->transcripts($entry, $alt_allele);
        my $consequence  = $transcript->{consequence};
        my @types        = split /\&/, $consequence;

        if (Genome::VariantReporting::Suite::Vep::SpliceNonsynonymousList::is_splice_site(@types)
            or Genome::VariantReporting::Suite::Vep::SpliceNonsynonymousList::is_non_synonymous(@types)) {
            $return_values{$alt_allele} = 1;
        }
        else {
            $return_values{$alt_allele} = 0;
        }
    }

    return %return_values;
}


sub vcf_id {
    return 'SPLICENONSYNONYMOUS';
}


sub vcf_description {
    return 'Variant hits either splice site or non-synonymous coding region';
}

1;




