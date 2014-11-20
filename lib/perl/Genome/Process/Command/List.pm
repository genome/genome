package Genome::Process::Command::List;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Process::Command::List {
    is => 'Genome::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1,
            value => 'Genome::Process'
        },
        show => { default_value => 'id,status,created_by,metadata_directory' },
    ],
    doc => 'list genome processes',
};


1;
