package Genome::VariantReporting::Reporter::AcmgReporter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Reporter::AcmgReporter {
    is  => ['Genome::VariantReporting::Reporter::FullReporter'],
    has => [
    ],
};

sub name {
    return 'acmg';
}

sub requires_interpreters {
    return (shift->SUPER::requires_interpreters, 'nhlbi');
}

sub headers {
    return (shift->SUPER::headers, 'All_MAF', 'AA_MAF', 'EU_MAF');
}

1;
