package Genome::Annotation::Expert::BamReadcount::Annotate;

use strict;
use warnings FATAL => 'all';

class Genome::Annotation::Expert::BamReadcount::Annotate {
    is => 'Genome::Annotation::Expert::CommandBase',
    has_input => [
        readcount_results => {
            is => 'Genome::Annotation::Expert::BamReadcount::RunResult',
            is_many => 1,
        },
    ],
};

sub execute {
    my $self = shift;

    $self->output_result(Genome::Annotation::Expert::BamReadcount::AnnotateResult->get_or_create($self->input_hash));
    return 1;
}

