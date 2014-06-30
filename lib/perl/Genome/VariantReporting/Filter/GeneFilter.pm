package Genome::VariantReporting::Filter::GeneFilter;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::VepConsequenceParser;

class Genome::VariantReporting::Filter::GeneFilter {
    is => ['Genome::VariantReporting::Component::Filter'],
    has => [
        gene_list => {
            is => 'String',
            doc => 'Name of the gene list to utilize for filtering',
            valid_values => ['acmg'],
        },
    ],
};

sub name {
    return 'gene';
}

sub requires_experts {
    return ('vep');
}

sub filter_entry {
    my $self = shift;
    my $entry = shift;

    my %return_values;
    my $vep_parser = new Genome::File::Vcf::VepConsequenceParser($entry->{header});
    my $desired_genes = $self->get_gene_list;

    for my $alt_allele ( @{$entry->{alternate_alleles}} ) {
        my ($transcript) = $vep_parser->transcripts($entry, $alt_allele);
        my $gene = $transcript->{'symbol'};
        if (defined $gene and exists $desired_genes->{uc($gene)}) {
            $return_values{$alt_allele} = 1;
        } else {
            $return_values{$alt_allele} = 0;
        }
    }

    return %return_values;
}

sub get_gene_list {
    my $self = shift;
    my $method = $self->gene_list . '_gene_list';

    if ($self->can($method)) {
        return $self->$method;
    } else {
        die $self->error_message('No method (%s) defined to get gene_list (%s)', $method, $self->gene_list);
    }
}

sub acmg_gene_list {
    my @genes = qw( ACTA2 ACTC1 APC APOB BRCA1 BRCA2 CACNA1S COL3A1 DSC2 DSG2 DSP FBN1 GLA KCNH2 KCNQ1
        LDLR LMNA MEN1 MLH1 MSH2 MSH6 MUTYH MYBPC3 MYH11 MYH7 MYL2 MYL3 MYLK NF2 PCSK9 PKP2 PMS2 PRKAG2
        PTEN RB1 RET RYR1 RYR2 SCN5A SDHAF2 SDHB SDHC SDHD SMAD3 STK11 TGFBR1 TGFBR2 TMEM43 TNNI3 TNNT2
        TP53 TPM1 TSC1 TSC2 VHL WT1
    );

    my %genes;
    map { $genes{$_} = 1 } @genes;
    return \%genes;
}

1;
