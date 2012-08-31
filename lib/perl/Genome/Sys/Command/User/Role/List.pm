package Genome::Sys::Command::User::Role::List;

use strict;
use warnings;
use Genome;

class Genome::Sys::Command::User::Role::List {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name => {
            is_constant => 1,
            value => 'Genome::Sys::User::Role',
        },
        show => { default_value => 'id,name' },
    ],
    doc => 'list user roles',
};

1;

