package Genome::Disk::Command::Assignment::List;

use strict;
use warnings;

use Genome;

class Genome::Disk::Command::Assignment::List {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1,
            value => 'Genome::Disk::Assignment',
        },
        show => { 
            default_value => 'disk_group_name,absolute_path,total_kb,unallocated_kb,percent_allocated' 
        },
    ],
};

1;
