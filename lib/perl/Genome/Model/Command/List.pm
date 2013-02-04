package Genome::Model::Command::List;
use strict;
use warnings;
use Genome;

class Genome::Model::Command::List {
    is => 'Genome::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::Model' 
        },
        show => { default_value => 'id,name,subject,processing_profile.' },
    ],
    doc => 'list genome models',
};

sub sub_command_sort_position { 3 }

1;

