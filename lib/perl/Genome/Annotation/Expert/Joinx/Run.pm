package Genome::Annotation::Expert::Joinx::Run;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Expert::Joinx::Run {
    is => 'Genome::Annotation::Expert::CommandBase',
    has_input => [
        known_variants => {
            is => 'Genome::Model::Build::ImportedVariationList',
            is_many => 1,
        },
        info_string => {
            is => 'Text',
        },
        version => {
            is => 'Text',
        },
    ],
};

sub name {
    'joinx';
}

sub execute {
    my $self = shift;

    $self->output_result(Genome::Annotation::Expert::Joinx::RunResult->get_or_create($self->input_hash));
    return 1;
}
