package Genome::VariantReporting::Reporter::VcfIndelsReporter;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::Writer;
use List::AllUtils qw(all);

class Genome::VariantReporting::Reporter::VcfIndelsReporter {
    is => 'Genome::VariantReporting::Reporter::VcfReporter',
};

sub name {
    return 'vcf-indels';
}

sub requires_interpreters {
    return qw(vcf-entry ft-keep contains-tag coverage-vaf genotype-vaf max-indel-size);
}

1;
