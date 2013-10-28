package Genome::Model::Tools::CopyCat::List;
use strict;
use warnings;
use Genome;

class Genome::Model::Tools::CopyCat::List {
    is => 'Genome::Object::Command::List',
    has => [
        subject_class_name => {
            is_constant => 1,
            value => 'Genome::Model::Tools::CopyCat::AnnotationData',
        },
    ],
    doc => 'List CopyCat annotation data sets',
};

1;
