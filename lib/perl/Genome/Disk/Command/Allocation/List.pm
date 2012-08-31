package Genome::Disk::Command::Allocation::List;

use strict;
use warnings;

use Genome;

class Genome::Disk::Command::Allocation::List {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::Disk::Allocation' 
        },
        show => { default_value => 'id,absolute_path,kilobytes_requested,owner_class_name,owner_id,creation_time,reallocation_time' },
    ],
    doc => 'lists Genome::Disk::Allocation objects',
};

1;

