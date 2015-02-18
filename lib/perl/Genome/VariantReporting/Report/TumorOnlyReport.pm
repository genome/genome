package Genome::VariantReporting::Report::TumorOnlyReport;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Report::TumorOnlyReport {
    is => 'Genome::VariantReporting::Report::WithHeader',
    doc => 'Basic variant information and annotation from vep, bam readcount, and population allele frequencies',
};

sub name {
    return 'tumor-only';
}

sub required_interpreters {
    return qw(position
    rsid
    caf
    vaf
    vep
    nhlbi
    1kg
    );
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
        sift
        polyphen
        condel
        dbSNP_caf
        NHLBI_All_MAF
        NHLBI_AA_MAF
        NHLBI_EU_MAF
        1kg-af
        vaf
        ref_count
        var_count
    /;
}

sub available_fields_dict {
    my $self = shift;
    my %dict = $self->SUPER::available_fields_dict;
    my $interpreter = $self->interpreters->{'vaf'};
    my @sample_names = $interpreter->sample_names;
    if (scalar(@sample_names) > 1) {
        die $self->error_message('Expecting only one sample to be used for the vaf interpreter. Got %s: %s', scalar(@sample_names), join(', ', @sample_names));
    }
    for my $field (qw(vaf ref_count var_count)) {
        $dict{$field} = {
            interpreter => 'vaf',
            field => $interpreter->create_sample_specific_field_name($field, $sample_names[0]),
        };
    }
    return %dict;
}

1;
