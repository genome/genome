package Genome::VariantReporting::Dbsnp::Run;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Dbsnp::Run {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Command',
    has_input => [
        dbsnp_build_id => {
            is => 'Text',
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
    'dbsnp';
}

sub result_class {
    'Genome::VariantReporting::Dbsnp::RunResult';
}
