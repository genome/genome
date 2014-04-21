package Genome::Annotation::Joinx::Run;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Joinx::Run {
    is => 'Genome::Annotation::CommandBase',
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

sub execute {
    my $self = shift;

    $self->output_result(Genome::Annotation::Joinx::RunResult->get_or_create($self->input_hash));
    return 1;
}
