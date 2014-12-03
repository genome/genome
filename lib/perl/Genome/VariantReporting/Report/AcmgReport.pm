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
    return (shift->SUPER::headers, 'All_MAF', 'AA_MAF', 'EU_MAF');
}

1;
