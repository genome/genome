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
        show => {
            value => 'id,reference_sequence,version,output_dir',
        }
    ],
    doc => 'List CopyCat annotation data sets',
};

sub _hint_string {
    return '';
}

1;
