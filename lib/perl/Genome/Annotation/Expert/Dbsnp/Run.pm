package Genome::Annotation::Expert::Dbsnp::Run;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Expert::Dbsnp::Run {
    is => 'Genome::Annotation::Expert::CommandBase',
    has_input => [
        known_variants => {
            is => 'Genome::Model::Build::ImportedVariationList',
            is_many => 1,
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

sub execute {
    my $self = shift;

    $self->output_result(Genome::Annotation::Expert::Dbsnp::RunResult->get_or_create($self->input_hash));
    return 1;
}
