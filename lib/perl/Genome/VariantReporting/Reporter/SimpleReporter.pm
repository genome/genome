package Genome::VariantReporting::Reporter::SimpleReporter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Reporter::SimpleReporter {
    is => 'Genome::VariantReporting::Reporter::WithHeader',
    has => [
    ],
};

sub name {
    return 'simple';
}

sub requires_interpreters {
    return qw(position vep ft contains-tag coverage vaf indel-size genotype);
}

sub headers {
    return qw/
        chromosome_name
        start
        stop
        reference
        variant
        transcript_name
        trv_type
        amino_acid_change
        default_gene_name
        ensembl_gene_id
        ft_string
        contains_tag
        coverage
        vaf
        genotype
        indel_size
    /;
}

1;
