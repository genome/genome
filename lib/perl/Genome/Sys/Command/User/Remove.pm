package Genome::Sys::Command::User::Remove;

use strict;
use warnings;
use Genome;
use Carp;

class Genome::Sys::Command::User::Remove {
    is => 'Command::V2',
    has => [
        users => {
            is => 'Genome::Sys::User',
            shell_args_position => 1,
            is_many => 1,
            doc => 'one or more users to be deleted',
        },
    ],
    doc => 'remove user(s) from the system',
};

sub _is_hidden_in_docs {
    return !Genome::Sys->current_user_is_admin;
}

sub execute {
    my $self = shift;

    my @names = map { $_->username } $self->users;

    for my $user ($self->users) {
        my $username = $user->username;
        my $rv = $user->delete;
        unless ($rv) {
            Carp::confess "Could not delete user $username";
        }
    }

    $self->status_message("Removed users: " . join(',', @names));
    return 1;
}

1;

