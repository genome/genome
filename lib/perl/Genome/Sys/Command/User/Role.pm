package Genome::Sys::Command::User::Role;

use strict;
use warnings;
use Genome;

class Genome::Sys::Command::User::Role {
    is => 'Command::Tree',
    doc => 'commands that work with user roles',
};

1;

