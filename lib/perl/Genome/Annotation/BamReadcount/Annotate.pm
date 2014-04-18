package Genome::Annotation::BamReadcount::Annotate;

use strict;
use warnings FATAL => 'all';

class Genome::Annotation::BamReadcount::Annotate {
    is => 'Genome::Annotation::Detail::Command',
    has_input => [
        readcount_results => {
            is => 'Genome::Annotation::BamReadcount::RunResult',
            is_many => 1,
        },
    ],
};

sub execute {
    my $self = shift;

    $self->output_result(Genome::Annotation::BamReadcount::AnnotateResult->get_or_create($self->input_hash));
    return 1;
}

