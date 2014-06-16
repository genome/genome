package Genome::VariantReporting::Expert::BamReadcount::Annotate;

use strict;
use warnings FATAL => 'all';

class Genome::VariantReporting::Expert::BamReadcount::Annotate {
    is => 'Genome::VariantReporting::Expert::CommandBase',
    has_input => [
        readcount_results => {
            is => 'Genome::VariantReporting::Expert::BamReadcount::RunResult',
            is_many => 1,
        },
    ],
};

sub result_class {
    'Genome::VariantReporting::Expert::BamReadcount::AnnotateResult';
}

1;
