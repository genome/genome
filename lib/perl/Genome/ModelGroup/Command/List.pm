package Genome::ModelGroup::Command::List;

use strict;
use warnings;

use Genome;

class Genome::ModelGroup::Command::List {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::ModelGroup' 
        },
        show => { default_value => 'id,name', doc => 'properties of the model-group to list (comma-delimited)', is_optional => 1 },
    ],
    doc => 'list model-groups',
};

1;

