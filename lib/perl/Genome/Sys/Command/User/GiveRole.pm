package Genome::Sys::Command::User::GiveRole;

use strict;
use warnings;

use Genome;
use Carp;

class Genome::Sys::Command::User::GiveRole {
    is => 'Command::V2',
    has => [
        role => {
            is => 'Genome::Sys::User::Role',
            shell_args_position => 1,
            doc => 'role to give to users',
        },
        users => {
            is => 'Genome::Sys::User',
            is_many => 1,
            shell_args_position => 2,
            doc => 'users to which the role will be given',
        },
    ],
    doc => 'provided user is given the provided role',
};

sub _is_hidden_in_docs {
    return !Genome::Sys->current_user_is_admin;
}

sub execute {
    my $self = shift;

    my $role = $self->role;
    for my $user ($self->users) {
        my $rv = $user->add_user_role($role);
        unless ($rv) {
            Carp::confess "Could not add role " . $role->name . " to user " . $user->username;
        }
    }

    $self->status_message("Added role " . $role->name . " to users: " . join(',', map { $_->username } $self->users));
    return 1;
}

1;

