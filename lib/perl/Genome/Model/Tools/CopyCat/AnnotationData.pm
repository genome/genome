package Genome::Model::Tools::CopyCat::AnnotationData;

use strict;
use warnings;

use Carp;
use Genome;


class Genome::Model::Tools::CopyCat::AnnotationData {
    is => 'Genome::SoftwareResult',

    has => [
        'reference_sequence' => {
            is => 'Genome::Model::Build::ReferenceSequence',
            is_input => 1,
        },

        'version' => {
            is => 'Integer',
            is_param => 1,
        },
    ],
};


sub annotation_data_path {
    my $self = shift;

    return $self->allocation->absolute_path;
}


1;
