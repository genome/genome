package Genome::VariantReporting::Suite::Vep::GeneFilter;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::VepConsequenceParser;
use Set::Scalar;

class Genome::VariantReporting::Suite::Vep::GeneFilter {
    is => ['Genome::VariantReporting::Framework::Component::Filter'],
    has => [
        gene_set => {
            is => 'String',
            doc => 'Name of the gene set to utilize for filtering',
            valid_values => ['acmg'],
        },
    ],
    doc => q{Filter out variants that aren't part of the specified gene set},
};

sub name {
    return 'gene';
}

sub requires_annotations {
    return ('vep');
}

sub filter_entry {
    my $self = shift;
    my $entry = shift;

    my %return_values;
    my $vep_parser = new Genome::File::Vcf::VepConsequenceParser($entry->{header});
    my $desired_genes = $self->get_gene_set;

    for my $alt_allele ( @{$entry->{alternate_alleles}} ) {
        my ($transcript) = $vep_parser->transcripts($entry, $alt_allele);
        my $gene = $transcript->{'symbol'};
        if ( defined $gene and $desired_genes->contains(uc $gene) ) {
            $return_values{$alt_allele} = 1;
        } else {
            $return_values{$alt_allele} = 0;
        }
    }

    return %return_values;
}

sub get_gene_set {
    my $self = shift;
    my $method = $self->gene_set . '_gene_set';

    if ($self->can($method)) {
        return Set::Scalar->new($self->$method);
    } else {
        die $self->error_message('No method (%s) defined to get gene_set (%s)', $method, $self->gene_set);
    }
}

# This list was obtained from http://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/
# GATA2, RUNX1, CEBPA added per dlarson
sub acmg_gene_set {
    return qw( ACTA2 ACTC1 APC APOB BRCA1 BRCA2 CACNA1S COL3A1 DSC2 DSG2 DSP FBN1 GLA KCNH2 KCNQ1
        LDLR LMNA MEN1 MLH1 MSH2 MSH6 MUTYH MYBPC3 MYH11 MYH7 MYL2 MYL3 MYLK NF2 PCSK9 PKP2 PMS2 PRKAG2
        PTEN RB1 RET RYR1 RYR2 SCN5A SDHAF2 SDHB SDHC SDHD SMAD3 STK11 TGFBR1 TGFBR2 TMEM43 TNNI3 TNNT2
        TP53 TPM1 TSC1 TSC2 VHL WT1 GATA2 RUNX1 CEBPA
    );
}

sub vcf_id {
    my $self = shift;
    return "GENE_IS_" . $self->gene_set;
}

sub vcf_description {
    my $self = shift;
    return "This transcript is part of the " . $self->gene_set . " gene set";
}

1;
