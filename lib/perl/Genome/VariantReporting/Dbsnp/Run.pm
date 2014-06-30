package Genome::VariantReporting::Dbsnp::Run;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Dbsnp::Run {
    is => 'Genome::VariantReporting::Component::Expert::Command',
    has_input => [
        known_variants => {
            is => 'Genome::Model::Build::ImportedVariationList',
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

sub execute {
    my $self = shift;

    $self->output_result(Genome::VariantReporting::Dbsnp::RunResult->get_or_create($self->input_hash));
    return 1;
}
