package Genome::Model::SomaticVariation::Command::List;

use strict;
use warnings;

use Genome;
use Command; 
use Data::Dumper;

class Genome::Model::SomaticVariation::Command::List {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::Model::SomaticVariation' 
        },
        show => { default_value => 'id,name,subject,tumor_model,normal_model' },
    ],
    doc => 'list somatic-variation genome models',
};

1;

