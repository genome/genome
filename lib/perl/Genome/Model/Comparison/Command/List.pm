package Genome::Model::Comparison::Command::List;

use strict;
use warnings;

use Genome;
use Command; 
use Data::Dumper;

class Genome::Model::Comparison::Command::List {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::Model::Comparison' 
        },
        show => { default_value => 'id,name,processing_profile,changes,from_models' },
    ],
    doc => 'list genome model comparisons',
};

1;

