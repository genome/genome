package Genome::Sys::Command::User::List;

use strict;
use warnings;
use Genome;
use Carp;

class Genome::Sys::Command::User::List {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name => {
            is_constant => 1,
            value => 'Genome::Sys::User',
        },
        show => { default_value => 'username,email,name' },
    ],
    doc => 'list users',
};

1;

