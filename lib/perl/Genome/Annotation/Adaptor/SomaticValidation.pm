package Genome::Annotation::Adaptor::SomaticValidation;

use strict;
use warnings;
use Genome;

class Genome::Annotation::Adaptor::SomaticValidation {
    is => 'Genome::Annotation::Adaptor',
    has_input => [
        build => {
            is => 'Genome::Model::Build::SomaticValidation',
        },
    ],
};

sub resolve_bam_results {
    my $self = shift;
    return [ $self->build->control_merged_alignment_result, $self->build->merged_alignment_result ];
}

1;
