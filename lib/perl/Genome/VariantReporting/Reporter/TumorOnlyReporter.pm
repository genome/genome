package Genome::VariantReporting::Reporter::TumorOnlyReporter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Reporter::TumorOnlyReporter {
    is => 'Genome::VariantReporting::Reporter::WithHeader',
};

sub name {
    return 'tumor-only';
}

sub requires_interpreters {
    return qw(position
    rsid
    gmaf
    vaf
    vep);
}

sub headers {
    return qw/
        chromosome_name
        start
        stop
        reference
        variant
        rsid
        transcript_name
        trv_type
        amino_acid_change
        default_gene_name
        gene_name_source
        ensembl_gene_id
        c_position
        canonical
        gmaf
        vaf
        ref_count
        var_count
    /;
}

1;
