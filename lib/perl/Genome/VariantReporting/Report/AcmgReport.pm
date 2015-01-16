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
    my @headers = grep {$_ ne 'caf' && $_ ne 'max_alt_af'} shift->SUPER::headers;
    push @headers, ('caf', 'max_alt_af', 'All_MAF', 'AA_MAF', 'EU_MAF');
    return @headers;
}

1;
