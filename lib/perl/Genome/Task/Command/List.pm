package Genome::Task::Command::List;

use warnings;
use strict;
use Genome;
use JSON::XS;
use File::Path;

class Genome::Task::Command::List {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::Task' 
        },
        show => { default_value => 'id,user_id,status,time_submitted', doc => 'properties of the model-group to list (comma-delimited)', is_optional => 1 },
    ],
    doc => 'list model-groups',
};

