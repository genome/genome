package Genome::VariantReporting::Reporter::VcfHomopolymerReporter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Reporter::VcfHomopolymerReporter {
    is => 'Genome::VariantReporting::Reporter::VcfReporter',
};

sub name {
    return 'vcf-homopolymer';
}

sub requires_interpreters {
    return qw(vcf-entry ft-keep contains-tag homopolymer);
}

1;
