package Genome::VariantReporting::Reporter::AcmgReporter;

use strict;
use warnings;
use Genome;
use List::Util qw( min );

class Genome::VariantReporting::Reporter::AcmgReporter {
    is  => ['Genome::VariantReporting::Reporter::FullReporter'],
    has => [
    ],
};

sub name {
    return 'acmg';
}

sub requires_interpreters {
    return qw(position vep info-tags variant-type min-coverage min-coverage-observed max-vaf-observed variant-callers many-samples-vaf rsid caf nhlbi);
}

sub headers {
    my $self = shift;
    my @headers = qw/
        chromosome_name
        start
        stop
        reference
        variant
        variant_type
        transcript_name
        trv_type
        trv_type_category
        amino_acid_change
        default_gene_name
        ensembl_gene_id
        inSegDup
        AML_RMG
        rsid
        caf
        All_MAF
        AA_MAF
        EU_MAF
        max_alt_af
        onTarget
        MeetsMinDepthCutoff
    /;

    push @headers, $self->_single_vaf_headers;

    push @headers, qw/
        min_coverage_observed
        max_normal_vaf_observed
        max_tumor_vaf_observed
        variant_callers
        variant_caller_count
    /;

    push @headers, $self->_per_library_vaf_headers;

    return @headers;
}

1;
