package Genome::Sys::Command::User::RemoveRole;

use strict;
use warnings;
use Genome;
use Carp;

class Genome::Sys::Command::User::RemoveRole {
    is => 'Command::V2',
    has => [
        role => {
            is => 'Genome::Sys::User::Role',
            shell_args_position => 1,
            doc => 'role to remove from provided users',
        },
        users => {
            is => 'Genome::Sys::User',
            shell_args_position => 2,
            is_many => 1,
            doc => 'users from which to remove the role',
        },
    ],
    doc => 'removes the provided role from the user',
};

sub _is_hidden_in_docs {
    return !Genome::Sys->current_user_is_admin;
}

sub execute {
    my $self = shift;

    for my $user ($self->users) {
        my $rv = $user->remove_user_role($self->role);
        unless ($rv) {
            Carp::confess "Could not remove role " . $self->role->name . " from user " . $self->user->username;
        }
    }

    $self->status_message("Successfully removed role " . $self->role->name . " from users: " .
        join(', ', map { $_->username } $self->users));
    return 1;
}

1;

