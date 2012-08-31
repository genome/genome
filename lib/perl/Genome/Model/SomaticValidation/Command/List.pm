package Genome::Model::SomaticValidation::Command::List;

use strict;
use warnings;

use Genome;
use Command; 
use Data::Dumper;

class Genome::Model::SomaticValidation::Command::List {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::Model::SomaticValidation' 
        },
        show => { default_value => 'id,name,subject,tumor_sample,normal_sample' },
    ],
    doc => 'list clinseq genome models',
};

1;

