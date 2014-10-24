package Genome::VariantReporting::Suite::Joinx::Run;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Suite::Joinx::Run {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Command',
    has_input => [
        vcf => {
            is => 'Path',
        },
        info_string => {
            is => 'Text',
        },
        joinx_version => {
            is => 'Text',
        },
    ],
};

sub name {
    return 'joinx';
}

sub result_class {
    'Genome::VariantReporting::Suite::Joinx::RunResult';
}
