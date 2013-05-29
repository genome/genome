package Genome::Model::ReferenceAlignment::Command::List;

use strict;
use warnings;

use Genome;
use Command; 
use Data::Dumper;

class Genome::Model::ReferenceAlignment::Command::List {
    is => 'Genome::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::Model::ReferenceAlignment' 
        },
        show => { default_value => 'id,name,subject,reference_sequence_build,processing_profile' },
    ],
    doc => 'list reference alignment genome models',
};

sub sub_command_sort_position { 3 }

1;

