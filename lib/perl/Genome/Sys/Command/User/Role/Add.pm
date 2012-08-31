package Genome::Sys::Command::User::Role::Add;

use strict;
use warnings;
use Genome;
use Carp;

class Genome::Sys::Command::User::Role::Add {
    is => 'Command::V2',
    has => [
        name => {
            is => 'Text',
            shell_args_position => 1,
            doc => 'name of new role',
        },
    ],
    doc => 'creates a new user role',
};

sub _is_hidden_in_docs {
    return !Genome::Sys->current_user_is_admin;
}

sub execute {
    my $self = shift;

    my $role = Genome::Sys::User::Role->create(
        name => $self->name,
    );
    unless ($role) {
        Carp::confess "Could not create role " . $self->name;
    }

    $self->status_message("Successfully created role " . $self->name);
    return 1;
}

1;

