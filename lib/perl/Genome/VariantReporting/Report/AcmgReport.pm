package Genome::VariantReporting::Report::AcmgReport;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Report::AcmgReport {
    is  => ['Genome::VariantReporting::Report::FullReport'],
    doc => 'Full report with additional population allele frequencies from NHLBI',
};

sub name {
    return 'acmg';
}

sub required_interpreters {
    return (shift->SUPER::required_interpreters, 'nhlbi');
}

sub headers {
    my @headers = grep {$_ ne 'dbSNP_caf' && $_ ne 'dbSNP_max_alt_af' && $_ ne 'max_tumor_vaf_observed'} shift->SUPER::headers;
    push @headers, ('dbSNP_caf', 'dbSNP_max_alt_af', 'NHLBI_All_MAF', 'NHLBI_AA_MAF', 'NHLBI_EU_MAF');
    return @headers;
}

1;
