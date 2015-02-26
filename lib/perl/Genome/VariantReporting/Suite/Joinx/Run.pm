package Genome::VariantReporting::Suite::Joinx::Run;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Suite::Joinx::Run {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Command',
    is_abstract => 1,
    has_planned_transient => [
        vcf => {
            is => 'Path',
            is_translated => 1,
            doc => 'Vcf File containing annotation',
        },
        info_string => {
            is => 'Text',
            doc => 'Field ids to embed from the annotation VCF. Use colons to separate multiple field descriptors.',
        },
        joinx_version => {
            is => 'Version',
            doc => 'joinx version to be used.',
        },
    ],
};

sub result_class {
    'Genome::VariantReporting::Suite::Joinx::RunResult';
}
