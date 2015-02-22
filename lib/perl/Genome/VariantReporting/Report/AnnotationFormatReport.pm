package Genome::VariantReporting::Report::AnnotationFormatReport;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Report::AnnotationFormatReport {
    is => 'Genome::VariantReporting::Report::WithHeader',
    doc => 'Report to mimic the legacy annotation format',
};

sub name {
    return 'annotation-format';
}

sub required_interpreters {
    return qw(
        position
        vep
    );
}

sub headers {
    return qw/
        chromosome_name
        start
        stop
        reference
        variant
        type
        gene_name
        transcript_name
        transcript_species
        transcript_source
        transcript_version
        strand
        transcript_status
        trv_type
        c_position
        amino_acid_change
        ucsc_cons
        domain
        all_domains
        deletion_substructures
        transcript_error
        default_gene_name
        gene_name_source
        ensembl_gene_id
        /;
}

sub unavailable_headers {
    return qw(
        type
        gene_name
        transcript_species
        transcript_source
        transcript_version
        strand
        transcript_status
        ucsc_cons
        domain
        all_domains
        deletion_substructures
        transcript_error
    );
}

1;
