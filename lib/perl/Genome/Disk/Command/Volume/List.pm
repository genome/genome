package Genome::Disk::Command::Volume::List;

use strict;
use warnings;

use Genome;

class Genome::Disk::Command::Volume::List {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1,
            value => 'Genome::Disk::Volume',
        },
        show => { 
            default_value => 'mount_path,disk_group_names,total_kb,percent_used,percent_allocated', 
        },
    ],
    doc => 'Lists Genome::Disk::Volume objects',
};

1;

